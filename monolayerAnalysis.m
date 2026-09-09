%% read in the data
tmp = readtable('thompsonLoringWaterFilms.xlsx');
ts = tmp{:,1:2:11};
fo = tmp{:,2:2:12};
mlVals = [170, 120, 53, 21, 5.8, 4.7];
rateTimes = [10, 30];

%% plot the experimental data
cs = parula(length(mlVals)); % insert your favourite colours here
ind = ts > rateTimes(1) & ts < rateTimes(2);
ps = gobjects(length(mlVals), 1);
figure;
for i = 1:length(mlVals)
    plot(ts(:,i), fo(:,i), 'LineWidth', 5, 'Color',...
                    cs(i,:), 'LineStyle', '-.');
    hold on;
    ps(i) = plot(ts(ind(:,i), i), fo(ind(:,i), i), 'LineWidth', 5, 'Color', cs(i,:));
end
l = legend(ps, num2str(mlVals'));
l.Location = 'northwest';
l.Title.String = 'Water film [ML]';
l.FontSize = 18;
xlim([0, 120]);
xlabel('Time [hr]');
ylabel('Carbonation [%]');

%% linearise the carbonation rate
rateStep = fo;
rateStep(~ind) = nan;
maxCarb = max(rateStep, [], 'omitnan');
minCarb = min(rateStep, [], 'omitnan');
totCarb = maxCarb - minCarb;
totCarb = totCarb / abs(diff(rateTimes)) / 3600 / 100; %hr to s, % to frac

%% query x values from 0.01 ML to 400 ML
qX = linspace(-4, 2.6021, 1000);
qX = 10.^qX;

% saturation function (Michaelis-Menten)
% using exponentials of b to guarantee positive values for b1 and b2
modelfun = @(b, x) exp(b(1)) * (x.^b(3)) ./ (exp(b(2)) + x .^ b(3));

% initial guesses
b1 = [-8, 2.5, -1];
f2 = fitnlm(table(mlVals', totCarb'), modelfun, b1);
[coffs, R, J, covB] = nlinfit(mlVals', totCarb', modelfun, b1);
% coffs and covB are used in the rate calculations function

%% generate random coeffs based on the covariance
n = 50000;
mcPars = mvnrnd(coffs, covB, n);

%% monte carlo fits – 1 fit for each set of coefficients. 
mcLines = nan(length(qX), n);
for i = 1:n
    mcLines(:,i) = feval(modelfun, mcPars(i,:), qX);
end

%% best fit line - define maximum rate at 400 ML 
bf = feval(modelfun, coffs, qX);
maxRate = bf(end);

%% quantiles of the lines
qnts = quantile(mcLines', [0.025, 0.16, 0.5, 0.84, 0.975]);

%% plot many lines and the data and the inset
cs = [0.3859, 0.4858, 0.7640];
figure;
hold on;
% % plot an example subset of the monte carlo lines
% for i = 1:ceil(n * 0.02)
%     plot(qX, mcLines(:,i) / maxRate, 'Color', [cs, 0.05]);
% end
scatter(mlVals, totCarb / maxRate, 100, [0 0 0], 'filled');
plot(qX, bf / maxRate, 'LineWidth', 3, 'LineStyle', '-', 'Color', [0 0 0]);
% plot(qX, qnts(1,:) / maxRate, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0]);
% plot(qX, qnts(5,:) / maxRate, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0]);
ylim([0, 1]);
xlim([0.01, 200]);
ylabel('Relative carbonation rate');
xlabel('Water film [ML]');

% inset plot with log scale
axes('Position', [0.45, 0.16, 0.4, 0.35]);
box on;
hold on;
plot(qX, bf / maxRate, 'LineWidth', 3, 'LineStyle', '-', 'Color', [0 0 0]);
% plot(qX, qnts(1,:) / maxRate, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0]);
% plot(qX, qnts(5,:) / maxRate, 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0 0]);
scatter(mlVals, totCarb / maxRate, 100, [0 0 0], 'filled');
xlim([0.9, 205]);
ax = gca;
set(ax, 'XScale', 'log');
grid off;
set(ax, 'LineWidth', 5);
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.TickLength = [0.03 0.0250];
ax.FontSize = 18;

