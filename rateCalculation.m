%% function for extrapolating our carbonation rates, dependent upon 
% temperature, partial pressure of CO2, water film thickness (in
% monolayers). mcRate is an optional return of n monte carlo simulations
% used to generate uncertainty estimates. 
function [carbRate, mcRate] = rateCalculation(temperature, pCO2, ml, n)
    ml(isnan(ml)) = 0;
    if nargin < 4
        n = 5000;
    end
    
    % import the fitting parameters 
    tFun = load('tempParameters.mat');
    mlFun = load('mlParameters.mat');
    B = mlFun.coffs;
    covB = mlFun.covB;
    mlModel = mlFun.modelfun;
    
    % generate n sets of the B parameters
    mcB = mvnrnd(B, covB, n);
    
    % generate n samples for each ml 
    relRates = nan(length(ml), n);
    for i = 1:n
        relRates(:,i) = feval(mlModel, mcB(i,:), ml);
    end
    
    % best-fit ML dependence
    bfML = feval(mlModel, B, ml);
    
    % max rate calculated as the best-fit parameters evaluated at 400 ML
    % thickness. 
    maxRate = feval(mlModel, B, 400);
    relRates = relRates / maxRate;
    bfML = bfML / maxRate;
    
    % now extract the temperature functions
    p = tFun.p;
    covP = tFun.covP;
    
    % generate n sets of pol parameters
    mcP = mvnrnd(p, covP, n);
    
    carbRates = nan(length(temperature), n);
    for i = 1:n
        carbRates(:,i) = polyval(mcP(i,:), 1 ./ temperature);
    end
    carbRates = exp(carbRates);
    
    % best-fit temperature relationship
    bfT = polyval(p, 1 ./ temperature);
    bfT = exp(bfT);
    
    fact = (sqrt(pCO2) / sqrt(150)) /... % pCO2
                                26.667 / 2; % 75µm -> 2 mm, fe(II) 0.4
    
    % mcRate is the monte-carlo simulations, assuming the rate-ML and rate-T 
    % are indepndent 
    mcRate = carbRates .* relRates * fact;
    
    % carbRate is the best fit
    carbRate = bfT.*bfML * fact;
end

