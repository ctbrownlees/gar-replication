%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk main replication script          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script replicates tables 4-6 of the paper.
%
%% Run setup 

global HOME;
run('code/setup');

% GaR Replication Analysis Options

% Beginning of out of sample period
os = 0.25;
% (The forecasting exercise starts at observation floor(os*T) where 
% T is the sample size)

% Vector of forecast horizons
H = [1 2 3 4]; 

% Nominal coverage level of GaR forecasts
coverage = 0.95;  

%% Run Analysis

% Create data
create_data(os);

% Run QR out-of-sample
qr_oos(H,coverage);

% Run GARCH out-of-sample
garch_oos(H,coverage); 

% Create tables for out-of-sample results
create_tables(H,coverage);
% Tables will be created in .csv files in HOME/tables. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
