%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Backtesting Global Growth-at-Risk Replication Script               % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script replicates the tables 1 to 6 reported on the paper.
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
    
addpath('code/matlab/utilities/')
addpath('code/matlab/insample/')
addpath('code/matlab/out-of-sample/')

% This block for LASSO estimation of quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also                                                                %
addpath('code/matlab/utilities/SDPT3-4.0/')                               %
addpath('code/matlab/utilities/SDPT3-4.0/Solver')                         %
addpath('code/matlab/utilities/SDPT3-4.0/Solver/Mexfun')                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Useful functions
scale_var = @(x) (x - mean(x))./sqrt(var(x));
%% Data Cleaning
% We leave 5% of the sample for ``training'' forecast errors from
% historical benchmarks.
is = 0.20;
os = 0.25;
% Create data
data_status = create_data(is,os);

H  = [1 2 3 4];
coverage = 0.95;
%% In-sample Results
% Table 1: QR Univariate
% Univariate models
models = {'nfci','ts','gf','cs','hp','sv','epu','cr','cg','wui','gpr'};
% Load some enviroment variables
ENV    = load(sprintf('%s/data/mfiles/data.mat',HOME));

for m = 1:length(models)
    try
        if strcmp(string(models(m)), 'gf')
            uni_status(m) = qr_univariate([] , [] , H , coverage , string( models( m )) );
        else
        var = readtable(sprintf('%s/data/clean/%s.csv', HOME , string( models( m ) )));
        try
            var_dates = datenum(var{:,2});
        catch % If date is in different formats, try this
            var_dates = datenum( cell2mat( var{:,2} ),...
    'YYYY-QQ');
        end
        var = scale_var( var{:,3:end} );
        var = var(ismember(var_dates, ENV.dates_nfci),:);
        var_dates = intersect(ENV.dates_nfci,var_dates);
        uni_status(m) = qr_univariate(var , var_dates , H , coverage , string( models( m )) );
        end
    catch
        sprintf('issue at %s',string(models(m)))
    end
end

% Table 2: QR Multivariate
qr_multi_status  = qr_multivariate( H , coverage );

% Table 3: GARCH
garch_is_status  = garch_insample( 0.95 );

%% Out-of-Sample Results
% QR out-of-sample
qr_os_status    = qr_oos(H,coverage,1);

% GARCH out-of-sample
garch_os_status = garch_oos( H , coverage); 

% Create tables for out-of-sample results
tables = create_tables(H , coverage) ;

