function create_data(is,os)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates the data.mat file that will be used throughout the code
% 
% INPUTS:
%   is   - Percentage of observations to be used for insample estimation of 
%           the historical BJPR.
%   os   - The percentage of observations to be used as insample for the
%           remaining models.
% 
% OUTPUTS:
%  This script creates a data.mat file in the data/inputs folder. 
%
% COMMENTS: 
% The BJPR is based on bootstrapping forecast errors. To generate
% forecast errors for the historical benchmark, we start forecasting
% historical BJPR at 'is'. Note also that is has to be smaller than os.
% All percentages are relative to the number of NFCI observations
% available.
%
% Authors: Christian Brownlees and Andre B.M. Souza
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global HOME % uses the global path defined in the setup file.

%% Load Data
data  = readtable( sprintf('%s/data/input/gdp.csv' , HOME ),'TreatAsEmpty' ,...
    { '.' , 'NA' } );
nfci_data  = readtable( sprintf('%s/data/input/nfci.csv', HOME), ...
    'TreatAsEmpty',{ '.' , 'NA' } );

dates_garch = data{ : , strcmp( data.Properties.VariableNames , 'Dates' ) };
dates_nfci = nfci_data{ : , 2 };

nfci       = nfci_data{ : , 3:size( nfci_data , 2 ) };

[common_dates, indexes_nfci, indexes_garch ] = intersect( dates_nfci,...
    dates_garch );
first_idx = find( sum( isnan( data{ : , 3:size( data , 2 ) } ), 2 ) == 0 ,...
    1 , 'first'); % Collect the first row for which there is no nan obs.

gdp_garch  = data{ first_idx : indexes_garch( end ), 3:size( data , 2 ) }; 
gdp_qr     = data{ ismember( data.Dates , dates_nfci ), 3:size( data , 2 ) };

dates_garch = datenum( cell2mat( dates_garch( first_idx:indexes_garch( end ) ) ),...
    'YYYY-QQ');
dates_nfci = datenum( cell2mat( dates_nfci ),...
    'YYYY-QQ');

%% Setting up parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T_nfci , N] = size( gdp_qr ); % We set sizes relative to NFCI
isnfci       = floor( is * T_nfci ); 
is           = sum( dates_garch <= dates_nfci( isnfci ) );
os           = floor( os * T_nfci ) ;
os_nfci      = os;
os_garch     = sum( dates_garch <= dates_nfci( os ) );
T            = size(gdp_garch,1);
%%  Exporting datafile   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(sprintf('%s/data/input/data.mat' , HOME ))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
