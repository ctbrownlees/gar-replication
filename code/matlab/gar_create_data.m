%% Create Data
% This script creates the required m-files and sets parameters.
%% Clearing everything
clear all
clc

%% Load Data
data  = readtable('../../data/clean/gdp.csv','TreatAsEmpty',{'.','NA'});
dates_garch = data{:,strcmp(data.Properties.VariableNames,'Dates')};
nfci_data  = readtable('../../data/clean/nfci.csv','TreatAsEmpty',{'.','NA'});
gdp_qr = data{end-175:end,3:size(data,2)};
nfci = nfci_data{:,3:size(nfci_data,2)};
dates       = dates_garch(end-175:end);
gdp_garch   = data{:,3:size(data,2)};

names = data.Properties.VariableNames(3:end);

% Fix sizes
[T , N]   = size(gdp_qr);

% Set balanced panel 

NAN = isnan(gdp_qr);
for i = 1:N
    lastNaN(i) = 0;
    if any(find(NAN(:,i),1,'last'))
        lastNaN(i) = find(NAN(:,i),1,'last');
    end
end
maxnan = 1;
for l = 1:T
    hasnan(l) = mean(any(isnan(gdp_garch(l:end,lastNaN<maxnan))))>0;
end

gdp_garch   = gdp_garch(sum(hasnan)+1:end,lastNaN<maxnan); % 
dates_garch = datenum(cell2mat(dates_garch(57:end)),'YYYY-QQ');
dates = datenum(cell2mat(dates(maxnan:T)),'YYYY-QQ');
gdp   = gdp_garch;
[T , N]   = size(gdp);

imf = {'KOR','AUS','BRA','CAN','CHL','CHI','FRA','DEU','IND','IDN','ITA',...
    'JPN','MEX','RUS','SAF','ESP','SWE','SWI','TUR','USA','GBR'};
imf_names = names(ismember(names,imf));
imf_sel   = find(ismember(names,imf));

%% Setting up parameters
is = floor(0.2*T) ; osFcast  =  floor(0.25*T) ;
isnfci = is ; osFcastnfci = osFcast;
is = sum(dates_garch <= dates(is)); osFcast = sum(dates_garch <= dates(osFcast));
Tn = T; T = size(gdp_garch,1);

%% Saving Names Table and data Mfile
%save('../data/mfiles/data.mat')
