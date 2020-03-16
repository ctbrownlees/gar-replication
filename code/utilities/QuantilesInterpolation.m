function [lc,sc,sh,df] = QuantilesInterpolation(qqTarg,QQ,lc0,sc0,sh0)
%% QuantilesInterpolation: Fit skewed t-distribution to target quantiles
% 
% Description: QuantilesInterpolation computes parameters for the skewed
% t-distribution to match a provided set of target quantiles. The
% parameters are chosen to minimize the squared distance between the target
% quantiles and the quantiles of the fitted skewed t-distribution.
% Optionally, initial conditions can be provided for the location, scale,
% and shape parameters.
% 
% The degrees of freedom are restricted to be integer-valued. To deal with
% this constraint, we search along a grid of values from 1 to 30 and
% optimize over the three continuous-valued parameters at each gridpoint.
% 
% Input arguments:
% - qqTarg : Vector with length = length(QQ) containing the target
%            quantiles to which the skewed t-distribution will be matched.
%            Only the 5th, 25th, 75th, and 95th quantiles will be used.
% - QQ : Vector of numbers between 0 and 1 (exclusive) indicating the
%        quantiles to which the entries of qqTarg correspond. This should
%        contain 0.05, 0.25, 0.50, 0.75, and 0.95.
% - lc0 : Initial condition for fitting the location parameter.
% - sc0 : Initial condition for fitting the scale parameter.
% - sh0 : Initial condition for fitting the shape parameter.
% 
% Output arguments:
% - lc : Location parameter for the fitted distribution.
% - sc : Scale parameter for the fitted distribution.
% - sh : Shape parameter for the fitted distribution.
% - df : Degrees of freedom parameter for the fitted distribution.
%% Set bounds and options for optimization
%     location  scale  shape
LB = [     -20,     0,   -30];
UB = [      20,    50,    30];

options =  optimset('MaxFunEvals', 1e+4,...
                    'MaxIter', 1e+4,...
                    'TolX', 1e-6,...
                    'TolCon', 1e-6,...
                    'Display', 'off',...
                    'LargeScale', 'on');
%% Locate target quantiles
[~, jq50] = min(abs(QQ - 0.50));
[~, jq25] = min(abs(QQ - 0.25));
[~, jq75] = min(abs(QQ - 0.75));
[~, jq05] = min(abs(QQ - 0.05));
[~, jq95] = min(abs(QQ - 0.95));

%% Set initial conditions for optimization (if not provided)
if nargin < 3
    iqn = norminv(0.75) - norminv(0.25);    
    lc0 = qqTarg(jq50);
    sc0 = (qqTarg(jq75) - qqTarg(jq25)) / iqn;
    sh0 = 0;
end

% Initial conditions
X0 = [lc0, sc0, sh0];

% Indices of target quantiles
Select = [jq05, jq25, jq75, jq95];

% Estimate approximation error for each possible value of the degrees of
% freedom parameter, optimizing over the other three continuous-valued
% parameters.
par = NaN(30, 3);
ssq = NaN(30, 1);
for df = 1:30
    [par(df, :), ssq(df)] = lsqnonlin(@(x) qqTarg(Select) - qskt(QQ(Select), x(1), x(2), x(3), df),...
                                      X0(1:3), LB, UB, options);
end
% Find degree of freedom value that provides the best fit, along with the
% three other parameters.
X = NaN(1, 4);
[~, X(4)] = min(ssq);
X(1:3) = par(X(4), :);

lc = X(1);
sc = X(2);
sh = X(3);
df = X(4);
