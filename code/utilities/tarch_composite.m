function [parameters,ht,zt,htfcast] = tarch_composite(data,shocks,p,o,q,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes GARCH models using composite likelihood. 
%
% INPUTS:
%   data   - T by N matrix of obvservations 
%   shocks - T by N vector of shocks; used to drive the Direct GARCH
%            for iterated GARCH, this is the same as data.
%   p      - order of the ARCH 
%   o      - order of the asymmetry; 0 is usual GARCH.
%   q      - order of the GARCH 
%   h      - forecast horizon
%
% OUTPUTS:
%   parameters  - Vector of estimated parameters 
%   ht          - Matrix of conditional volatilities 
%   zt          - Matrix of standardized residuals 
%   htfcast     - Forecast of the conditional variances 
%   
% COMMENTS:
%  This function estimates the parameters of the following equation 
%  
%  h_{i,t} = var * (1 - alpha - beta) + alpha .* data^2 + beta h_{i,t-1}
%
%  where var is the variance of the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options
options = optimset('fminunc');
options = optimset(options,'TolFun',1e-005);
options = optimset(options,'TolX',1e-005);
options = optimset(options,'MaxFunEvals',10000);
options = optimset(options,'Display','Off');
% Do a grid search over parameters
if o>=1
    par = [ [0.01:0.005:0.105]' , [0.01:0.005:0.105]' ,[0.95:-0.01:0.76]'];
    for j=1:size(par,1)
        for k=1:size(par,1)
            LLG(j,k) = tarch_composite_LL([0.01 par(j,2),par(k,3)],data,shocks,p,o,q,h);
        end
    end
else
    par = [ [0.01:0.01:0.2]' ,[0.96:-0.01:0.77]'];
    for j=1:size(par,1)
        for k=1:size(par,1)
            LLG(j,k) = tarch_composite_LL([par(j,1),par(k,2)],data,shocks,p,o,q,h);
        end
    end
end

[mr , mc ] = find(LLG==min(LLG(:)));
% First Minimization on the min of grid
if o>=1
   [parameters,LL] =  fminunc('tarch_composite_LL',[0.01 par(mr),par(mc)],options,data,shocks,p,o,q,h);
else
   [parameters,LL] =  fminunc('tarch_composite_LL',[par(mr,1),par(mc,2)],options,data,shocks,p,o,q,h);
end

% Robustness against typical parameters
if o>=1
   [parameterscheck,LLcheck] =  fminunc('tarch_composite_LL',[0.01 0.15,0.9],options,data,shocks,p,o,q,h);
else
   [parameterscheck,LLcheck] =  fminunc('tarch_composite_LL',[0.05,0.90],options,data,shocks,p,o,q,h);
end
if LL > LLcheck
    parameters = parameterscheck;
end
% Set final output
[~ , ht,zt,htfcast] = tarch_composite_LL(parameters,data,shocks,p,o,q,h);
