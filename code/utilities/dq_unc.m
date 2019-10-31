function[pval] = dq_unc(p,h,nahead)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the p-value of a ``direct'' DQ test based on
%  a constant. 
%
% USAGE:
%   [pval] = dq_hits( p , h , nahead)
%
% INPUTS:
%   p        - Quantile of interest 
%   h        - Hit sequence
%   nahead   - Forecast horizon; used to adjust the covariance of hits. 
%
% OUTPUTS:
%   pval     - p-value of the null hypothesis that b_i = 0 for all i
%              in the following equation
%
%   h_{t+nahead} = b_0 + e_{t+nahead} 
%   
% COMMENTS:
% We need to think about the covariance estimation under M dependence.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n    = length(h);
  y    = h(1:n)-p;
  X    = ones(size(y));
  k    = size(X,2);
  if nahead > 1
    [b , covmat] = olsnw_dq(y,X,0,nahead-1,p);
    cinv = inv(covmat);
    stat = b'*cinv*b ;
  else
    b    = (X'*X)\(X'*y);
    k    = size(X,2);
    stat = ( b' * X' * X * b )./(p*(1-p));
  end
  pval = 1-chi2cdf(stat, k) ;
end
