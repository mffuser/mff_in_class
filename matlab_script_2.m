%% download DAX data

tickSymb = '^GDAXI';

dateBeg = '01011990';
dateEnd = '01072015';
dateEnd = datestr(today, 'ddmmyyyy');

%%

daxCrude = hist_stock_data(dateBeg, dateEnd, tickSymb);

%%

nObs = length(daxCrude.Date);
serialDates = zeros(nObs, 1);

tic;
for ii=1:nObs
    serialDates(ii, 1) = datenum(daxCrude.Date{ii, 1}, 'yyyy-mm-dd');
end
time = toc;

%%

tic;
serialDates = datenum(daxCrude.Date, 'yyyy-mm-dd');
time2 = toc;

%%

dax.dates = flipud(serialDates);
dax.prices = flipud(daxCrude.AdjClose);

%%

plot(dax.prices)

%%

figure('position', [50 50 1200 600])
subplot(1, 2, 1);
plot(dax.prices)

subplot(1, 2, 2);
plot(dax.dates, dax.prices)
datetick 'x'
xlabel('dates')
ylabel('prices')
title('historic DAX values')

set(gca, 'xLim', [dax.dates(1) dax.dates(end)])


%% get maximum

[val, ind] = max(dax.prices);

%%

maxVal = max(dax.prices);
isFulFilled = zeros(nObs, 1);

for ii=1:nObs
   if dax.prices(ii, 1) == maxVal
       isFulFilled(ii, 1) = 1;
   end
end

%%
maxInd = find(isFulFilled);
dax.prices(maxInd)
dax.dates(maxInd)

%%

plot(dax.dates, dax.prices)
datetick 'x'
xlabel('dates')
ylabel('prices')
title('historic DAX values')

set(gca, 'xLim', [dax.dates(1) dax.dates(end)])

hold on;
plot(dax.dates(maxInd), dax.prices(maxInd), '.r')

%% logical indexing

matr = [1 2 3 4 5 6];

% check if entries are greater than 3
greaterThan = matr > 3;

matr(greaterThan)

%% 

%equalThree = matr == 3;
matr(matr ~= 3)

%%

matr(matr == 3) = 8

%%


log1 = rand(3) > 0.5
log2 = rand(3) > 0.5

AND = (log1 & log2)
OR = (log1 | log2)
NOTBOTH = ~AND

%%

a = [0 3 1 0 2 9];
b = 3;
assert(isequal(nearZero(a),b));


%% volatility clusters

dax.disRet = 100*(dax.prices(2:end) - dax.prices(1:(end-1))) ...
    ./ dax.prices(1:(end-1));

dax.retDates = dax.dates(2:end);

ax(1) = subplot(2, 1, 1);
plot(dax.retDates, dax.prices(2:end))
datetick 'x'
set(gca, 'xLim', [dax.retDates(1) dax.retDates(end)])

ax(2) = subplot(2, 1, 2);
plot(dax.retDates, dax.disRet)
datetick 'x'
set(gca, 'xLim', [dax.retDates(1) dax.retDates(end)])

linkaxes([ax(1) ax(2)], 'x')

%% 

hist(dax.disRet, 30)

meanRet = mean(dax.disRet);
stdDev = std(dax.disRet);

yLimits = get(gca, 'yLim');
line([meanRet meanRet], yLimits, 'Color', 'r')
line((meanRet+2*stdDev)*[1, 1], yLimits, 'Color', 'r')
line((meanRet-2*stdDev)*[1, 1], yLimits, 'Color', 'r')
text(meanRet+2*stdDev, yLimits(end)/2, '2 standard dev')

shg

%% NaNs

a = [1; NaN; 3]
b = [1; NaN; 3]

a == b

%%

NaN > NaN
3 + NaN


%%

a = [12; 10; NaN; 14];

a(2:end) - a(1:end-1)

%% anonymous function

plusThree = @(x)x+3;

%% download data

tickSymb1 = '^GDAXI';
tickSymb2 = '^GSPC';

dateBeg = '01011990';
dateEnd = '01072015';
dataOut = hist_stock_data(dateBeg, dateEnd, ...
    tickSymb1, tickSymb2);

%% 

tickSymbs = {'^GDAXI', '^GSPC'};

dateBeg = '01011990';
dateEnd = '01072015';

data = getPrices(dateBeg, dateEnd, tickSymbs);

%%

dateBeg = '01011990';
dateEnd = '01072011';

% download data of all components: dax_comp is structure array
daxComp = {'ADS.DE', 'ALV.DE',...
    'BAS.DE', 'BAYN.DE', 'BEI.DE', 'BMW.DE', 'CBK.DE', 'DAI.DE', ...
    'DB1.DE',...
    'DBK.DE', 'DPW.DE', 'DTE.DE', 'EOAN.DE', 'FME.DE', 'FRE.DE',...
    'HEI.DE', 'HEN3.DE', 'IFX.DE', 'LHA.DE', 'LIN.DE', 'MAN.DE',...
    'MEO.DE', 'MRK.DE', 'MUV2.DE', 'RWE.DE', 'SAP', 'SDF.DE',...
    'SIE.DE', 'TKA.DE', 'VOW3.DE', '^GDAXI'};

profile on

data = getPrices(dateBeg, dateEnd, daxComp);





