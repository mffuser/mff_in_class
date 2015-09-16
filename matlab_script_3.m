%% Monte Carlo Simulation: uniform distribution
% 
% % init params
% n = 1000;
% nBins = 20;
% 
% simU = rand(n, 1);
% 
% hist(simU, nBins)
% 
% expNum = n/nBins;
% 
% line([0 1],[expNum expNum], 'Color', 'r')
% 
% %% Monte Carlo: normal distribution
% n2 = 10;
% 
% mu = 4;
% sigma = 1;
% 
% scatter(zeros(n2, 1), simU(1:n2), '.r')
% hold on
% 
% grid = 0:0.1:8;
% normVals = normcdf(grid, mu, sigma);
% plot(grid, normVals)
% 
% X = norminv(simU, mu, sigma);
% 
% for ii=1:n2
%    line([0 X(ii)], [simU(ii) simU(ii)], 'Color', 'r',...
%        'LineStyle', ':')
%    line([X(ii) X(ii)], [0 simU(ii)], 'Color', 'r',...
%        'LineStyle', ':')
% end
% 
% 
% %%
% 
% figure('position', [50 50 1200 600])
% 
% subplot(1, 2, 1)
% [counterU locationsU] = hist(simU, nBins);
% 
% width = diff(locationsU(1:2));
% 
% bar(locationsU, counterU / (n*width), 1)
% 
% subplot(1, 2, 2)
% [counterX locationsX] = hist(X, nBins);
% 
% widthX = diff(locationsX(1:2));
% bar(locationsX, counterX / (n*widthX), 1)
% 
% hold on;
% grid = mu-4:0.1:mu+4;
% plot(grid, normpdf(grid, mu, sigma), '-r')
% 
% %%
% 
% n = 1000;
% mu = 2;
% sigma = 1;
% 
% points = normrnd(mu, sigma, n, 1);
% 
% lh1 = (sigma^2*2*pi)^(-0.5)*exp((-0.5)*(points-mu).^2/sigma^2);
% 
% lhs = normpdf(points, mu, sigma);
% 
% llh1 = log(normpdf(points, mu, sigma));
% llh2 = -0.5*log(sigma^2*2*pi) - 0.5*(points-mu).^2/sigma^2;
% 
% llh = sum(llh2);
% 
% %% anonymous function to calculate llh
% nllh = @(params, data)sum(0.5*log(params(2)^2*2*pi)+...
%     0.5*(data - params(1)).^2/params(2)^2);
% 
% 
% %% optimization
% 
% params0 = [0 1];
% 
% lb = [-inf 0.1];
% ub = [inf inf];
% 
% %opt = optimset('algorithm', 'interior-point'
% 
% tic;
% paramsHat = fmincon(@(params)nllh(params,points),...
%     params0, [], [], [], [], lb, ub);
% time1 = toc;
% 
% paramsHat2 = fmincon(nllh,...
%     params0, [], [], [], [], lb, ub,[],[],points);
% 
% tic;
% paramsHat3 = fmincon(@(params)normlike(params, points),...
%     params0, [],[],[],[],lb,ub);
% time2 = toc;
% 
% paramsHat4 = fmincon(@normlike,...
%     params0, [], [], [], [], lb, ub,[],[],points);
% 
% %%
% 
% n = 1000;
% nu = 4;
% 
% X = tinv(rand(n,1), nu);
% 
% lhs = @(NU, Points)gamma((NU+1)/2)*(NU*pi)^(-0.5)/gamma(NU/2)*...
%     (1+Points.^2/NU).^(-(NU+1)/2);
% 
% tpdf(X(1),4)
% lhs(4, X(1))
% 
% nllh = @(nu, x)-sum(log(lhs(nu, x)));
% 
% % init guess
% param0 = 8;
% 
% lb = 0.1;
% ub = inf;
% 
% nuHat = fmincon(nllh, param0, [],[],[],[],lb,ub,...
%     [],[],X);
% 
% %%
% 
% dateBeg = '01011990';
% dateEnd = datestr(today, 'ddmmyyyy');
% tickSymb = {'^GDAXI'};
% 
% % get prices
% daxPrices = getPrices(dateBeg, dateEnd, tickSymb);
% 
% % calculate returns
% logRetsTable = price2retWithHolidays(daxPrices);
% 
% % extract returns as vector (not table)
% logRets = 100*logRetsTable{:, :};
% dats = datenum(logRetsTable.Properties.RowNames, 'yyyy-mm-dd');
% 
% plot(dats, logRets)
% datetick 'x'
% 
% % create anonymous likelihood function
% lhs = @(params, x)normpdf(x, params(1), params(2));
% 
% % create anonymous negative log likelihood function
% nllh = @(params, x)-sum(log(lhs(params,x)));
% nllh = @(params, x)-sum(log(normpdf(x, params(1), params(2))));
% % nllh = @(params)-sum(log(normpdf(logRets, params(1), params(2))));
% 
% % set optimization parameters (initial value, parameter bounds,...)
% params0 = [0 1];
% 
% lb = [-inf 0.0001];
% ub = [inf inf];
% 
% % conduct optimization
% paramsHat = fmincon(@(params)nllh(params, logRets),params0,...
%     [],[],[],[],lb,ub);
% 
% [muHat, sigmaHat] = normfit(logRets)
% 
% %% QQ-plot
% 
% nObs = length(logRets);
% sortedLogRets = sort(logRets);
% 
% alphas = (1:nObs)/(nObs+1);
% 
% normQuants = norminv(alphas, muHat, sigmaHat);
% 
% minVal = floor(min([logRets(:); normQuants(:)]));
% maxVal = ceil(max([logRets(:); normQuants(:)]));
% 
% plot(sortedLogRets, normQuants, '.')
% set(gca, 'xLim', [minVal maxVal], 'yLim', [minVal maxVal])
% axis square
% hold on;
% 
% line([minVal maxVal],[minVal maxVal])
% 
% %%
% 
% [counterRet locationsRet] = hist(logRets, 30);
% width = locationsRet(2) - locationsRet(1);
% 
% %% 
% 
% grid = floor(min(logRets(:))):0.1:ceil(max(logRets(:)));
% yVals = normpdf(grid, muHat, sigmaHat);
% 
% bar(locationsRet, counterRet/(nObs*width))
% hold on;
% plot(grid, yVals, 'r-')
% 
% set(gca, 'xLim', [grid(1) -3])
% 
% 
% 
% %% VaR: normal distribution
% 
% alphaLevels = 1 - [0.95, 0.99, 0.995];
% 
% var_norm = norminv(alphaLevels, muHat, sigmaHat);
% 
% %% Var: t distribution
% 
% nllh = @(nu, data)-sum(log(tpdf(data, nu)));
% 
% param0 = 4;
% lb = 0.5;
% ub = inf;
% 
% nuHat = fmincon(@(nu)nllh(nu, logRets), param0, [], [],...
%     [], [], lb, ub);
% 
% var_t = tinv(alphaLevels, nuHat);
% 
% %% Var: historical simulation
% 
% var_hs = quantile(logRets, alphaLevels);
% 
% %% all var values
% 
% [var_norm; var_t; var_hs]
% 
% %% backtesting
% 
% exceedFrequ = @(var_Vals)...
%     sum(repmat(logRets,1,3) < repmat(var_Vals, nObs, 1), 1)/nObs;
% 
% [alphaLevels;
%     exceedFrequ(var_norm);
%     exceedFrequ(var_t);
%     exceedFrequ(var_hs)]
% 
% %% Plot exceedances
% 
% var_Val = var_t(2);
% 
% figure('position', [50 50 1200 600])
% 
% scatter(dats, logRets, '.')
% datetick 'x'
% set(gca, 'xLim', [dats(1) dats(end)])
% 
% hold on;
% 
% 
% line([dats(1) dats(end)], var_Val*[1 1], 'Color', 'black')
% 
% exceedances = logRets < var_Val;
% dats_exceed = dats(exceedances);
% logRets_exceed = logRets(exceedances);
% 
% plot(dats_exceed, logRets_exceed, '.r')
% 
% %% autocorrelation
% 
% nLags = 20;
% 
% autoCorrCoef = zeros(1, nLags);
% for ii=1:nLags
%     autoCorrCoef(1, ii) = ...
%         corr(logRets(1:end-ii), logRets(1+ii:end));
% end
% 
% stem(1:nLags, autoCorrCoef, '.r', 'MarkerSize', 12)
% set(gca, 'yLim', [-0.5 1])
% 
% %%
% autocorr(logRets)
% 
% %% Autoregressive process
% 
% a = 0.8;
% n = 10000;
% sigma = 0.8;
% y0 = 0;
% 
% 
% epsilons = normrnd(0, sigma, n, 1);
% y = zeros(n, 1);
% y(1, 1) = y0;
% 
% for ii=2:n
%     y(ii, 1) = a*y(ii-1, 1) + epsilons(ii);
% end
% 
% figure('position', [50 50 1200 600])
% subplot(1, 2, 1)
% plot(y(1:150))
% 
% subplot(1, 2, 2)
% autocorr(y)
% 
% %%
% 
% a1 = 0;
% a2 = 0;
% a3 = -0.05;
% a4 = 0.05;
% 
% y0 = 0;
% y1 = 0;
% y2 = 0;
% y3 = 0;
% 
% epsilons = normrnd(0, sigma, n, 1);
% y = zeros(n, 1);
% y(1:4, 1) = [y0; y1; y2; y3];
% 
% for ii=5:n
%     y(ii) = a1*y(ii-1) + a2*y(ii-2) + a3*y(ii-3) + a4*y(ii-4) + epsilons(ii);
% end
% 
% figure('position', [50 50 1200 600])
% subplot(1, 2, 1)
% autocorr(y)
% 
% subplot(1, 2, 2)
% autocorr(logRets)
% 
% 
% %%
% 
% figure('position', [50 50 1200 600])
% subplot(1, 2, 1)
% autocorr(y.^2)
% 
% subplot(1, 2, 2)
% autocorr(logRets.^2)
% 

%% local variances

windowSize = 300;
empStds = zeros(nObs - windowSize + 1, 1);

for ii=windowSize:nObs
    empStds(ii-windowSize +1) = std(logRets((ii-windowSize+1):ii)); 
end

figure('position', [50 50 1200 600])
subplot(1,2,1)
plot(logRets)
subplot(1,2,2)
plot(empStds)


%% GARCH

% initial 
y0 = 0;
sigma0 = 1;

sampSize = 10000;

% garch parameters
garch = 0.9;
arch = 0.05;
k = 0.05;

% preallocation
y = zeros(sampSize, 1);
sigmas = zeros(sampSize, 1);

epsilons = randn(sampSize, 1);

y(1) = y0;
sigmas(1) = sigma0;

for ii=2:sampSize
    sigmas(ii) = ...
        sqrt(k + garch*sigmas(ii-1).^2 + arch*y(ii-1).^2);
    
    y(ii) = epsilons(ii)*sigmas(ii);
end


figure('position', [50 50 1200 800])

subplot(3, 2, 1:2)
plot(y(1:600))

subplot(3, 2, 3:4)
plot(sigmas(1:600))

subplot(3,2,5)
autocorr(y)

subplot(3,2,6)
autocorr(y.^2)

        
%% retrieve sigma values

% initial guess
sigmaHat0 = 2;

retrieveSigmas = zeros(numel(y), 1);
retrieveSigmas(1) = sigmaHat0;

for ii=2:numel(y)
    retrieveSigmas(ii) = sqrt(k + garch*retrieveSigmas(ii-1).^2 +...
        arch*y(ii-1).^2);
end



%% plotting

plot(sigmas)
hold on;
plot(retrieveSigmas, '-r')


%% calculate log likelihood

nllh = sum(0.5*log(retrieveSigmas.^2*2*pi) + ...
    0.5*(y.^2./retrieveSigmas.^2));


nllh2 = garchNllh([k garch arch], y)


%% fit garch

params0 = [0.1 0.4 0.4];

lb = [0.01 0.01 0.01];
ub = [1 1 1];

A = [0 1 1];
b = 1;

opt = optimset('algorithm', 'sqp');

paramsHat = fmincon(@(params)garchNllh(params, logRets), params0,...
    A, b, [], [], lb, ub, [],  opt);


%% retrieve real sigmas
sigmaHat0 = 1;

retrieveSigmas = zeros(numel(logRets), 1);
retrieveSigmas(1) = sigmaHat0;

for ii=2:numel(logRets)
    retrieveSigmas(ii) = sqrt(k + garch*retrieveSigmas(ii-1).^2 +...
        arch*logRets(ii-1).^2);
end

%%

figure('position', [50 50 1200 600])
subplot(2,1,1)
plot(dats, logRets)
datetick 'x'

subplot(2,1,2)
plot(dats, retrieveSigmas)
datetick 'x'































































