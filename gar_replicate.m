%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk main replication script          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the tables reported on the paper.
% It is organized in sections, each of which runs a function that saves
% a .mat file and outputs a status flag.
%% Clearing
clear all; clc
%% Global Variable definitions and paths
global HOME
dr = dir('.');
dr = {dr.name};
if ~sum(contains(dr,'code'))==1
    cd ../../ % Move to landing folder, if not there
    HOME = cd ;
else
    HOME = cd;
end
addpath('code/')    
addpath('code/utilities/')

% This block for LASSO estimation of quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also                                                                %
addpath('code/utilities/SDPT3-4.0/')                               %
addpath('code/utilities/SDPT3-4.0/Solver')                         %
addpath('code/utilities/SDPT3-4.0/Solver/Mexfun')                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GaR Replication Exercise Options
% We leave 5\% of the sample for ``training'' forecast errors from
% historical benchmarks. I DON'T UNDERSTAND WHAT THIS MEANS

is       = 0.20;      % NOT SURE WHAT THIS IS
os       = 0.25;      % NOT SURE WHAT THIS IS
H        = [1 2 3 4]; % Horizons
coverage = 0.95;      % Nominal coverage level of GaR forecasts

%% Run

% Create data
create_data(is,os);

% Run QR out-of-sample
qr_oos(H,coverage,1);

% Run GARCH out-of-sample
garch_oos( H , coverage); 

% Create tables for out-of-sample results
create_tables(H , coverage) ;

