%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk main replication script          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates tables  4-6 of the paper.
%% Run setup script
global HOME;
run('code/setup');

% GaR Replication Analysis Options
is       = 0.20;      % See create_data.m for more information. 
os       = 0.25;      % Out-of-sample forecasting starts at 25\% of the
% available observations for QR.
H        = [1 2 3 4]; % Vector of forecast horizons
coverage = 0.95;      % Nominal coverage level of GaR forecasts

%% Run Analysis

% Create data
create_data(is,os);

% Run QR out-of-sample
qr_oos(H,coverage,1);

% Run GARCH out-of-sample
garch_oos( H , coverage); 

% Create tables for out-of-sample results
create_tables(4, coverage) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
