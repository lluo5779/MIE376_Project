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
yearlyReturn.Properties.RowNames = cellstr(dates((size(dates,1)-size(yearlyReturn,1)+1):end, :));

monthlyReturn = readtable('monthly_ret.xlsx'); %xlsread('Project1_Data_adjClose_yearly_ret','Project1_Data_adjClose_yearly_r');
%monthlyReturn = monthlyReturn(:,2:end);
monthlyReturn.Properties.RowNames = cellstr(dates(1:size(monthlyReturn,1), :));

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

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Determine hyperparameters based on adjClose
[NoTotalDates, NoAssets] = size(adjClose);

% yearly training returns and prices from 2012 Jan to 2014 Dec
trainingReturns_yr = table2array( yearlyReturn( calStart <= dates & dates <= calEnd, :) );
lastYearTrainingReturns_yr = table2array( yearlyReturn( calEnd-calmonths(12)+days(1) <= dates & dates <= calEnd, :) );

% yearly testing returns from 2015 Jan to 2015 Dec
%testingReturns = table2array( yearlyReturn( testStart <= dates & dates <= testEnd, :) );

% weekly testing returns from 2015 Jan to 2015 Dec
testingReturns_weekly = table2array( returns( testStart <= dates & dates <= testEnd, :) );

% monthly training returns and prices from 2012 Jan to 2014 Dec
trainingReturns_mth = table2array( monthlyReturn( calStart <= dates & dates <= calEnd, :) );
lastYearTrainingReturns_mth = table2array( monthlyReturn( calEnd-calmonths(12)+days(1) <= dates & dates <= calEnd, :) );

% monthly testing returns from 2015 Jan to 2015 Dec
testingReturns = table2array( monthlyReturn( calStart <= dates & dates <= calEnd, :) );

% Data Mean Year
mu_entire_training_yr = geomean(trainingReturns_yr, 1) - 1;
mu_last_year_yr = geomean(lastYearTrainingReturns_yr, 1) - 1;

% Data Mean Month
mu_entire_training_mth = geomean(1+trainingReturns_mth, 1) - 1;
mu_last_year_mth = geomean(1+lastYearTrainingReturns_mth, 1) - 1;

% Data Variance Year
cov_entire_training_yr = cov(trainingReturns_yr-1);
cov_last_Year_yr = cov(lastYearTrainingReturns_yr-1);

% Data Variance Month
cov_entire_training_mth = cov(trainingReturns_mth);
cov_last_Year_mth = cov(lastYearTrainingReturns_mth);

no_scenarios = 10;
generated_ret = mvnrnd(mu_entire_training_mth, cov_entire_training_mth, no_scenarios); % 20 assets = 20 rows; scenarios in columns
exp_generated_ret = geomean(1+generated_ret,1)-1; % for VSS first term

%% CALL FUNCTION for Weights
[x,fval] = Solver(1+generated_ret, 1000, 1000)
fval = -1*fval;
[x_expected_scenarios, fval_expected_scenarios] = Solver(exp_generated_ret, 1000, 1200); %VSS first term


for i = 1:no_scenarios
    weights{i} = zeros(NoAssets,1);
    weights_exp_scenarios = zeros(NoAssets,1);
end
for i = 1:no_scenarios
    weights{i}(:,1) = x((NoAssets*i+1):NoAssets*(i+1), 1);%[x(1:NoAssets, 1); x((NoAssets*i+1):NoAssets*(i+1), 1)];
end

weights_exp_scenarios = x((NoAssets+1):NoAssets*2, 1);

%% 
for s = 1:no_scenarios
    portfRet(:,s) = testingReturns_weekly*weights{s};
end

expected_portfRet_across_scenarios = mean(portfRet);
portfRet_avg_scenarios = testingReturns_weekly*weights_expected_scenarios;

portfValue = weights*testingReturns_weekly';

plot(portfValue);

% y1[2]- y2[2]- y3[2]- y1[2]+ y2[2]+ y3[2]+ y1[3]- y2[3]- y3[3]- y1[3]+ y2[3]+ y3[3]+ 
% Aeq = [eye(3) eye(3) -1*eye(3) eye(3) zeros(3, no_var-12); % s = 1 and j = 1 2 3
%     eye(3) eye(3) zeros(3,6) -1*eye(3) eye(3) zeros(3, no_var-18); % s = 2 and j = 1 2 3
%     eye(3) eye(3) zeros(3,12) -1*eye(3) eye(3)] % s = 3 and j = 1 2 3







