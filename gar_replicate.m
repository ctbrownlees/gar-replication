%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk main replication script          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script replicates tables 4-6 of the paper.
%
%% Run setup 

global HOME; % Define the home path for the script as the landing folder.
run('code/setup');

% GaR Replication Analysis Options

% Beginning of out of sample period
os = 0.25;
% (The forecasting exercise starts at observation floor(os*T) where 
% T is the sample size)

% Vector of forecast horizons
H = [1 2 3 4 ]; 

% Nominal coverage level of GaR forecasts
coverage = [0.95 0.75 0.25 0.05];  
%% Run Analysis

% Create data
create_data(os);

% Run QR out-of-sample
qr_oos(H,coverage);

% Run Panel QR out-of-sample
panel_qr_oos(H,coverage);

% Run GARCH out-of-sample
garch_oos(H,coverage); 

% Constructing Distributions
distributions_oos(H);

% Create tables for out-of-sample results
create_tables(H,coverage);

% Tables will be created in .csv files in HOME/tables. 
create_tables_panel(H,coverage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
