function[pval] = aug_dq( p , h , x , nahead)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the p-value of a ``direct'' DQ test, where the
% information set againts which optimality is being tested includes
% the last 4 lags of x.
%
% USAGE:
%   [pval] = aug_dq( p , h , x , nahead)
%
% INPUTS:
%   p        - Quantile of interest 
%   h        - Hit sequence
%   x        - Variable included in the DQ test. 
%   nahead   - Forecast horizon; used to adjust the covariance of hits. 
%
% OUTPUTS:
%   pval     - p-value of the null hypothesis that b_i = 0 for all i
%              in the following equation
%
%   h_{t+nahead} = b_0 + b_1 x_{t-1} + ... + b_4x_{t-4} + e_{t+nahead} 
%   
% COMMENTS:
% We need to think about the covariance estimation under M dependence.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n    = length(h);
  y    = h(5:n)-p;
  X    = double([ones(size(y)), x(4:(n-1)) , x(3:(n-2)) , x(2:(n-3)) , x(1:(n-4))]);
  k    = size(X,2);  
 try
   if nahead > 1
    [b , covmat] = olsnw_dq(y,X,0,nahead-1,p);
    cinv   = inv(covmat);
    stat = b'* cinv * b;
  else
    b    = (X'*X)\(X'*y);
    k    = size(X,2);
    stat = ( b' * X' * X * b )./(p*(1-p));
  end
  pval = 1-chi2cdf(stat, k) ;
   catch
  pval = 0 ;
  end
end
 
