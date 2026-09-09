%% read in the data from O'Connor 2004
tmp = readtable('oconnor04appendix.xls', 'sheet', 'combined');

%% columns for time, final carbonation, temperature, co2, and grainsize
ts = tmp.Time_H_1;
rxTot = tmp.ExtentOfCarbonation__4;
temp = tmp.Temp__C_1 + 273;
co2 = tmp.PCO2_Atm;
sz = tmp.size_Um;

%% twin sisters olivine  
tsoInd = strcmp(tmp{:,2}, 'TSO');

%% index for twin sisters, 150 bar, and grainsize < 75µm
qInd = tsoInd & co2 == 150 & sz == -75;
x = temp(qInd);
y = rxTot(qInd) ./ ts(qInd);
y = y / 3600; % rate per second
y(y < 0) = 0;

% take only the low-temperature end (< 450 K) for arrhenius calculation
ltInd = qInd & temp < 450;
arrRate = rxTot(ltInd) ./ ts(ltInd);
arrRate = arrRate / 3600; % from per hr to per second.
arrTemp = temp(ltInd);
arrRate(arrRate < 0) = 0;

% query temperatures
qTemp = 180:2:550;

%% fit arrhenius - p and covP are used in the rate calculation function
ind = arrRate > 0;
[p, S] = polyfit(1 ./ arrTemp(ind), log(arrRate(ind)), 1);
pR = S.R;
covP = (inv(pR)*inv(pR)')*S.normr^2/S.df; % covariance matrix of p

%% generate random cofffs based on the covariance
n = 50000;
mcPars = mvnrnd(p, covP, n);
% p and covP are used in the rate calculation function. 

%% monte carlo fits – 1 fit for each set of coefficients. 
mcLines = nan(length(qTemp), n);
for i = 1:n
    mcLines(:,i) = polyval(mcPars(i,:), 1 ./ qTemp);
end
mcLines = exp(mcLines);

%% best fit line
bf = polyval(p, 1 ./ qTemp);
bf = exp(bf);

%% quantiles - 1 sig, 2 sig
qnts = quantile(mcLines', [0.025, 0.16, 0.5, 0.84, 0.975]);

%% plot the best fit line and the quantiles 
c1 = [81, 118, 100] / 255;

figure;
% LEFT TILE: 1/T x axis
tiledlayout(1, 2, 'TileSpacing', 'tight');

nexttile();
plot(1 ./ qTemp, bf, 'LineWidth', 3, 'Color', [0 0 0]); % convert to per second
hold on;
plot(1 ./ qTemp, qnts(1,:), 'LineWidth', 2, 'Color', [0.5 0.5 0.5],...
                                'LineStyle', '--');
plot(1 ./ qTemp, qnts(5,:), 'LineWidth', 2, 'Color', [0.5 0.5 0.5],...
                                'LineStyle', '--');
scatter(1 ./ x, y, 100, c1, 'filled',...
                         'LineWidth', 2, 'MarkerEdgeColor', [0 0 0]);

xlabel('1 / T [K^{-1}]');
ylabel('Reaction rate [1/s]');
set(gca, 'YScale', 'log');
xlim([1.86, 3] * 10^-3);
pbaspect([1, 1, 1])
grid off;

% RIGHT TILE: LINEAR T, extrapolation to low T
nexttile();
yticks([]);
yyaxis right;
plot(qTemp, bf, 'LineWidth', 3, 'Color', [0 0 0]);
hold on;
plot(qTemp, qnts(1,:), 'LineWidth', 2, 'Color', [0.5 0.5 0.5],...
                                'LineStyle', '--');
plot(qTemp, qnts(5,:), 'LineWidth', 2, 'Color', [0.5 0.5 0.5],...
                                'LineStyle', '--');
scatter(x, y, 100, c1, 'filled',...
                         'LineWidth', 2, 'MarkerEdgeColor', [0 0 0]);
set(gca, 'YScale', 'log');
xlim([180, 435]);
ylim([10^-12, 1.5 * 10^-1]);
xlabel('T [K]');
pbaspect([1, 1, 1])
grid off;
    