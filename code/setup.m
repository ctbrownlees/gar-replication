%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script sets up the required environment variables and paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clearing
clear all; clc
%% Global Variable definitions and paths
global HOME
dr = dir('.');
dr = {dr.name};
if ~sum(contains(dr,'code'))==1
    cd ../   % Move to landing folder, if not there
    HOME = cd ;
else
    HOME = cd;
end
addpath('code/')    
addpath('code/utilities/')

% This block for LASSO estimation of quantiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See also   https://faculty.fuqua.duke.edu/~abn5/belloni-software.html   %                              
addpath('code/utilities/SDPT3-4.0/')                                      %
addpath('code/utilities/SDPT3-4.0/Solver')                                %
addpath('code/utilities/SDPT3-4.0/Solver/Mexfun')                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
