function create_data(is,os)
global HOME
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

%% Setting up parameters
[T_nfci , N]   = size( gdp_qr ); % We set sizes relative to NFCI, the shortest panel
isnfci         = floor( is * T_nfci ); 
is             = sum( dates_garch <= dates_nfci( isnfci ) );
os             = floor( os * T_nfci ) ;
os_nfci        = os;
os_garch       = sum( dates_garch <= dates_nfci( os ) );
T              = size(gdp_garch,1);
%%  Exporting datafile
save(sprintf('%s/data/input/data.mat' , HOME ))
