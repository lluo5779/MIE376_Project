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
% monthlyReturn.Properties.RowNames = cellstr(dates(1:size(monthlyReturn,1), :));
avgMonthlyReturn = table2array(readtable('averaged_monthly_ret.xlsx'));

% Identify the tickers and the dates 
tickers = yearlyReturn.Properties.VariableNames';

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start of in-sample calibration period 
calStart = datetime('2013-01-01');
calEnd   = calStart + calmonths(12*2) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2015-01-01');
testEnd   = testStart + calmonths(12) - days(8);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Determine hyperparameters based on adjClose and monthlyReturn
[NoTotalDates, NoAssets] = size(adjClose);
[NoTotalMonths, NoAssets] = size(monthlyReturn);

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

no_splits = 50;
generated_ret = mvnrnd(mu_entire_training_mth, cov_entire_training_mth, no_splits); % 20 assets = 20 rows; scenarios in columns
exp_generated_ret = geomean(1+generated_ret,1)-1; % for VSS first term

b = 1000;
%G = [500, 900, 1000, 1100, 1200, 1500, 2000, 10000]
G = [950, 1000, 1050, 1100, 1200];

%% 3. Stochastic Program for varying parmeter
clear weights weights_mean weights_exp_scenarios fval_actual_scenarios fval_avg_scenarios


% Run solver for different G values
for i = 1:size(G, 2)
    % Call solver function to optimize weights for Stochastic Program
    [x,fval] = Solver_mat(1+generated_ret, b, G(1,i));
    fval_actual_scenarios(i,:) = -1*(fval);
    
    % Extract weights corresponding to each rollout of scenarios
    for j = 1:no_splits
        weights{i,j} = [x((1:NoAssets), 1) x((NoAssets*j+1):(NoAssets*(j+1)), 1)];
    end
    
    % Average of weights from stochastic problem
    weights_mat = cat(3, weights{i,:});
    weights_mean{i} = mean(weights_mat, 3);
    
    % Determine the weights for expected return = the true sample mean
    [x_expected_scenarios, fval_expected_scenarios] = Solver_mat(1+exp_generated_ret, b, G(1,i)); %VSS first term
    weights_avg_scenarios{i} = [x_expected_scenarios((1:NoAssets), 1) x_expected_scenarios(((NoAssets+1):NoAssets*2), 1)];
    fval_avg_scenarios(i,:) = (fval_expected_scenarios)*-1;
end



% 'weights{3,1}', 'weights{3,2}', 'weights{3,3}', 'weights{3,4}', 'weights{3,5}'
%% 4. Deterministic Program for t = 0 to t = 1

clear x_optimal_det det_value

% First we find the deterministic problem using expected returns from
% historical data

% Objective function with surplus reward of 1 and slack penalty of 4
f = [ones(1,NoAssets) -1 4];

A = [ones(1,NoAssets) 0 0];
b = 1000;

lb = zeros(1, NoAssets+2);

for i = 1:size(G,2)
    % Determine optimal x values for varying G
    beq = G(1,i);
    Aeq1 = [(avgMonthlyReturn((end-1), :)+1) -1 1];
    [x_optimal_det{i,1}, det_value{i,1}] = linprog(f, A, b, Aeq1, beq, lb, []);
    Aeq2 = [(avgMonthlyReturn(end, :)+1) -1 1];
    [x_optimal_det{i,2}, det_value{i,2}] = linprog(f, A, b, Aeq2, beq, lb, []);
end


save('very_important_table_det.mat', 'x_optimal_det', 'det_value')
%% 5. Value of Stochastic Solution

%weights_avg_scenarios
vss = fval_actual_scenarios-fval_avg_scenarios;
save('vss.m', 'fval_actual_scenarios','fval_avg_scenarios','vss')
%% 6. EVIP

% EVIP = 

%% 7. Portfolio Weights and Calculation of Portfolio Return for Out-of-Sample Periods
clear portfRet

for i = 1:size(G,2)
    fig1 = figure(i);

    X = weights_mean{i}(:,1)';
    labels = tickers;
    ax1 = subplot(1,2,1);
    p{(i-1)*2+1} = pie(ax1,X);
    title(ax1,strcat('G = ', num2str(G(1,i)),'at t = 0'),  'FontSize', 10);
    legend(labels, 'Location', 'bestoutside', 'Orientation', 'vertical');
    
    Y = weights_mean{i}(:,2)';
    ax2 = subplot(1,2,2);
    p{(i-1)*2+2} = pie(ax2,Y);
    title(ax2,strcat('G = ', num2str(G(1,i)),'at t = 1'), 'FontSize', 10);
    legend(labels, 'Location', 'bestoutside', 'Orientation', 'vertical');
    
    print(fig1,strcat('pi_', num2str(i)),'-dpng','-r0');
end


%% 8. Visualization of Portfolio Return for each Scenarios for Out-of-Sample Months
clear legend_text
% counter = 0
fig2 = figure(2)

mon_ret_for_testing = (1+avgMonthlyReturn((end-1):end,:))
for i = 1:size(G,2)
    
    for s = 1:no_splits
        if mod(s,2)==1
            t1_ret = mon_ret_for_testing(1,:)*weights{i,s}(:,1)/sum(weights{i,s}(:,1));
            t2_ret = mon_ret_for_testing(2,:)*weights{i,s}(:,2)/sum(weights{i,s}(:,2))*t1_ret;
            portfRet{i,s} = [t1_ret;t2_ret];
            plot([ones(1,1); portfRet{i,s}]'*1000)
            hold on
        end
%         counter = counter + 1
    end
end
% legend(legend_text)
title('Out-of-Sample: Portfolio Values for all Scenarios Rollouts', 'FontSize', 10);
ylabel('Portfolio Value','interpreter','latex','FontSize',12);
xlabel('Period')

print(fig2,'portfolio_value_all_scenarios','-dpng','-r0');
%% 9. Visualization of Portfolio Return for Averaged Scenarios for Out-of-Sample Months

fig3 = figure()
portfRet_mat = cat(3, portfRet{:});
expected_portfRet_across_scenarios = mean(portfRet_mat, 3);
weights_mean_mat = cat(3, weights_mean{:});
for i =1:size(G,2)
    t_ret = (1+avgMonthlyReturn((end-1):end,:))
    portfRet_testing(i,1) = (1+avgMonthlyReturn((end-1),:))*(weights_mean{i}(:,1)/sum(weights_mean{i}(:,1)));
    portfRet_testing(i,2) = portfRet_testing(i,1)*(1+avgMonthlyReturn((end),:))*(weights_mean{i}(:,2)/sum(weights_mean{i}(:,2)));
    portfRet_testing_det(i,1) = t_ret(1,:)*(x_optimal_det{i,1}(1:NoAssets)/sum(x_optimal_det{i}(1:NoAssets)));
    portfRet_testing_det(i,2) = portfRet_testing_det(i,1)*t_ret(2,:)*(x_optimal_det{i,2}(1:NoAssets)/sum(x_optimal_det{i}(1:NoAssets)));
    
end
% portfValue = weights*testingReturns_weekly';

plot([ones(size(G,2),1) portfRet_testing]'*1000);
hold on
plot([ones(size(G,2),1) portfRet_testing_det]'*1000, '--');

ylabel('Portfolio Value','interpreter','latex','FontSize',12);

print(fig3,'portfolio_value_avg_scenarios','-dpng','-r0');

% y1[2]- y2[2]- y3[2]- y1[2]+ y2[2]+ y3[2]+ y1[3]- y2[3]- y3[3]- y1[3]+ y2[3]+ y3[3]+ 
% Aeq = [eye(3) eye(3) -1*eye(3) eye(3) zeros(3, no_var-12); % s = 1 and j = 1 2 3
%     eye(3) eye(3) zeros(3,6) -1*eye(3) eye(3) zeros(3, no_var-18); % s = 2 and j = 1 2 3
%     eye(3) eye(3) zeros(3,12) -1*eye(3) eye(3)] % s = 3 and j = 1 2 3







