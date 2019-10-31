%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk main replication script          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the tables reported on the paper.
% It is organized in sections, each of which runs a function that saves
% a .mat file and outputs a status flag.
%% Run setup script
setup;

%% Data Cleaning
% We leave 5\% of the sample for ``training'' forecast errors from
% historical benchmarks.
is = 0.2 ;
os = 0.25;

% Create data.mat file that will be used throughout
create_data(is,os);

% Define horizons of interest and coverage levels
H  = [1 2 3 4];
coverage = 0.95;

%% Out-of-Sample Results
% QR out-of-sample
qr_oos(H,coverage,1);

% GARCH out-of-sample
garch_oos( H , coverage); 

% Create tables for out-of-sample results
create_tables(H , coverage) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
