%% read the input files in
t = readtable('dryPerplex_minModes.tab', 'FileType', 'text');
tComp = readtable('dryPerplex_totalComp.tab', 'FileType', 'text');
ss = readtable('dryPerplex_solidSolution.tab', 'FileType', 'text');

% rename relevant variable names 
tLabs = t.Properties.VariableNames;
tLabs = strrep(tLabs, 'f_CO2', 'fCO2');
tLabs = strrep(tLabs, 'f_O2', 'fO2');
tLabs = strrep(tLabs, 'O_HP_', 'ol');
tLabs = strrep(tLabs, 'Opx_HP_', 'opx');
tLabs = strrep(tLabs, 'iron', 'fe');
tLabs = strrep(tLabs, 'Cc_AE_', 'carb');
tLabs = strrep(tLabs, 'q', 'qz');
tLabs = strrep(tLabs, 'ta', 'tlc');
tLabs = strrep(tLabs, 'mt', 'mag');

%% normalise to the sum of MgO + FeO + SiO2 (rock weight = 1)
tCompOrig = sum(tComp{:,3:5}, 2);
tNorm = tComp{:,3:end};
tNorm = tNorm ./ tCompOrig;
% calculate the 'air' needed to create that composition. Assuming
% atmosphere of 10^-2 O2 and 10^-5 CO2. 
a2rock = tNorm(:,[6,5]) ./ [10.^-2, 10.^-5];
aNeeded = max(a2rock, [], 2);

%% define white background colours 
whites = ones(40, 3);

%% calculate the vector angle for streamline calculation
vectAngle = tNorm(:,[6,5]) ./ (10.^[tComp{:,1}, tComp{:,2}]);
vectAngle = -vectAngle ./ sum(abs(vectAngle), 2);
arrowScale = [0.2, 1.6];
normAngle = vectAngle .* arrowScale;

%% reshape all the composition points into a square, and the vectors of gas
% conspumption
m = sqrt(length(tComp{:,1}));
sqX = linspace(min(tComp{:,1}), max(tComp{:,1}), m);
sqY = linspace(min(tComp{:,2}), max(tComp{:,2}), m);
[sqX, sqY] = meshgrid(sqX, sqY);
sqU = reshape(normAngle(:,1), [m, m])';
sqV = reshape(normAngle(:,2), [m, m])';
sqAN = reshape(aNeeded, [m, m])'; % square aNeeded
sqMins = reshape(t{:,3:end}, [m, m, width(t)-2]); % mineral modes in square
for i = 1:length(sqMins(1,1,:))
    sqMins(:,:,i) = sqMins(:,:,i)';
end
sqSS = reshape(ss{:,3:end}, [m, m, width(ss)-2]); % solidSolution in square
for i = 1:length(sqSS(1,1,:))
    sqSS(:,:,i) = sqSS(:,:,i)';
end


%% calculate the streamline of the vector field - to get the reaction path
f = figure;
sl = streamline(sqX, sqY, sqU, sqV, -2, -5, [0.005, 10000000]);
slx = sl.XData;
sly = sl.YData;
slx = slx(~isnan(slx));
sly = sly(~isnan(sly));
close(f);

%% manually add the final corner of the olivine field (the streamline often 
% stops right before the end – comment this out to see the result)
olFin = [-8.356, -93.5];
slx = [slx, linspace(slx(end), olFin(1), 5)];
sly = [sly, linspace(sly(end), olFin(2), 5)];

%% INTERPOLATE VARIOUS METRIX ALONG THE STREAMLINE
%% take a moving min to smooth the air-rock ratio (i.e. remove the 
% excursions from the main reaction path). 
iAN = interp2(sqX, sqY, sqAN, slx', sly', 'nearest');
iAN = movmin(iAN, [500000, 0]);
% need a point that represents the unaltered phase (0 air:rock doesn't plot
% on the log x axis)
iAN(end) = 0.01;


%% interp the mineral phases along the reaction path
iMin = nan(length(slx), size(sqMins, 3));
for i = 1:size(sqMins, 3)
    tmp = interp2(sqX, sqY, sqMins(:,:,i), slx, sly, 'nearest');
    iMin(:,i) = tmp;
end

%% interp the solid-solution composition along the reaction path.
iSS = nan(length(slx), size(sqSS, 3)); 
for i = 1:size(sqSS, 3)
    tmp = interp2(sqX, sqY, sqSS(:,:,i), slx', sly', 'nearest');
    iSS(:,i) = tmp;
end


%% only plot values at unique air-rock ratios. 
[~, ind] = unique(iAN);

% plot the mineral abundance
figure;
tiledlayout(2, 1, 'TileSpacing', 'compact');
cs = parula(width(iMin));

nexttile;
minHands = gobjects(width(iMin), 1);
minBool = true(length(minHands),1);
for i=1:width(iMin)
    if all(isnan(iMin(ind,i)))
        minBool(i) = false;
    end
    hold on;
    minHands(i) = stairs(iAN(ind), iMin(ind, i) / 100, 'Color', cs(i,:), 'LineWidth', 5);
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([20, 4000]);
ylim([0.004, 1.2]);
xticklabels([]);
minLabs = tLabs(3:end);
legend(minHands(minBool), minLabs(minBool), 'Location', 'southwest');
ylabel('Mineral Fraction');
grid off;

% plot the solid solution
nexttile;
ssHands = gobjects(width(iSS)-1, 1);
cOrder = [1, 2, 6];
for i=1:width(iSS)
    hold on;
    ssHands(i) = stairs(iAN(ind), iSS(ind, i), 'Color',...
                                        cs(cOrder(i),:), 'LineWidth', 5);
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([20, 4000]);
ylim([0.48, 1.02]);
yticks(0.5:0.1:1);
grid off;
xlabel('air:rock');
ylabel('Mg/(Mg+Fe)');

%%
cs = parula(100); % choose your favourite colours here. 

% ind to plot only some subset of the vector field. You can make this
% denser or sparser by changing dens 
dens = 51;
ind = 1:dens:(m*m);

% plumx is a custom function for plotting perplex pseudosections. 
figure;
[~, ps, ts, f] = plumx(t, 'ExpansionValue', 0.6,...
                'SmoothValue', 0.2, 'LineWidth', 2.5,...
                 'MaxFontSize', 20, 'MinFontSize', 10,...
                 'TextOverflow', false, 'MineralLabels', tLabs,...
                 'PhaseColors', whites, 'LineColor', [0.2 0.2 0.2],...
                  'Numerate', true);
hold on;
colormap(cs);
lineCols = iAN;
lineCols(lineCols < 1) = nan; % omit the added point (0.01) and any zeros
scatter(slx, sly, 30, log10(lineCols), 'filled');
quiver(sqX(ind), sqY(ind), sqU(ind), sqV(ind), 'off',...
                        'ShowArrowHead', 'off', 'LineWidth', 1.5, 'Color',...
                            [0.7 0.7 0.7]);
                        
% creat the legend for the text that was enumerated.                        
tb = false(size(ts));
for i = 1:length(ts)
    if ~isnan(str2double(ts{i}.String))
        tb(i) = true;
    end
end    
legend(ps(tb), f(tb), 'Location', 'westoutside', 'FontSize', 16);
cb = colorbar;
cb.Label.String = 'log10(air:rock)';
xlabel('fCO2');
ylabel('fO2');

% these figures need to be finalised in illustrator, as one might expect.
% In particular, the vector field needs to have arrow heads added manually.
% ALWAYS COMPARE THE PHASE BOUNDARIES TO THE NATIVE PSSECT FROM PERPLE_X

