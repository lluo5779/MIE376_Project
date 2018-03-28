clc
clear all
format long

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project1_Data_adjClose.csv');
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
dates = datetime(adjClose.Date);
size_adjClose = size(adjClose);
adjClose = adjClose(:,2:size_adjClose(2));

factorRet = readtable('Project1_Data_FF_factors.csv'); %, 'ReadRowNames', true);
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));

yearlyReturn = readtable('yearly_ret.xlsx'); %xlsread('Project1_Data_adjClose_yearly_ret','Project1_Data_adjClose_yearly_r');
yearlyReturn = yearlyReturn(:,2:end);
yearlyReturn.Properties.RowNames = cellstr(dates((size(dates,1)-size(yearlyReturn,1)+1):end, :))

% Identify the tickers and the dates 
tickers = yearlyReturn.Properties.VariableNames';

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start of in-sample calibration period 
calStart = datetime('2013-01-01');
calEnd   = calStart + calmonths(12*2) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2015-01-01');
testEnd   = testStart + calmonths(12) - days(8);

% Initialize testing returns
testingReturns = table2array( returns( testStart <= dates & dates <= testEnd, :) );

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Determine hyperparameters based on adjClose
[NoTotalDates, NoAssets] = size(adjClose);

% period returns and prices from 2012 Jan to 2014 Dec
trainingReturns = table2array( yearlyReturn( calStart <= dates & dates <= calEnd, :) );
lastYearTrainingReturns = table2array( yearlyReturn( calEnd-calmonths(12)+days(1) <= dates & dates <= calEnd, :) );

% testing returns from 2015 Jan to 2015 Dec
testingReturns = table2array( yearlyReturn( calStart <= dates & dates <= calEnd, :) );

% Data Mean
mu_entire_training = geomean(1+trainingReturns, 1) - 1;
mu_last_year = geomean(1+lastYearTrainingReturns, 1) - 1;
                                 
% Data Variance
cov_entire_training = cov(trainingReturns);
cov_last_Year = cov(lastYearTrainingReturns);

no_scenarios = 3;
generated_ret = mvnrnd(mu_entire_training, cov_entire_training, no_scenarios);

%% CALL FUNCTION for Weights
weights = rand([1,20])

%%
portfRet = weights*testingReturns';








