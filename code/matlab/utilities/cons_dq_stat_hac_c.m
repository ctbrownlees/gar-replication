function[dq_t] = cons_dq_stat_hac(p,h,nahead)  
  n    = length(h);
  y    = h(1:n)-p;
  X    = ones(size(y));
  k    = size(X,2);
  if nahead > 1
    [b , ~ ,~ , covmat] = olsnw_dq(y,X,0,nahead-1,p);
    cinv   = inv(covmat);
    stat = b'*cinv*b ;
  else
    b    = (X'*X)\(X'*y);
    k    = size(X,2);
    stat = ( b' * X' * X * b )./(p*(1-p));
  end
  dq_t = 1-chi2cdf(stat, k) ;
end
