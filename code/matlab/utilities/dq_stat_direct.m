function[dq_t] = dq_stat_direct(p,h,nahead)  
  n    = length(h);
  y    = h(4+nahead:n)-p;
  X    = double([ones(size(y)), h(4:(n-nahead)) , h(3:(n-1-nahead)) , h(2:(n-2-nahead)) , h(1:(n-3-nahead))]);
  b    = (X'*X)\(X'*y);
  k    = size(X,2);
  stat = (b(2:k)' * X(:,2:k)' * X(:,2:k) * b(2:k) )./(p*(1-p));
%  dq_t = 1-chi2cdf(stat, k-1) ;
  dq_t = 1-chi2cdf(stat, k) ;
end
 
