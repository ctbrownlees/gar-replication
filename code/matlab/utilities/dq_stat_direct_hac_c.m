function[dq_t] = dq_stat_direct_hac_c(p,h,nahead)  
  n    = length(h);
  y    = h(4+nahead:n)-p;
  X    = double([ones(size(y)), h(4:(n-nahead)) , h(3:(n-1-nahead)) , h(2:(n-2-nahead)) , h(1:(n-3-nahead))]);
  k    = size(X,2);
  try
  if nahead > 1
    [b , ~ ,~ , covmat] = olsnw_dq(y,X,0,nahead-1,p);  
    cinv   = inv(covmat);
    stat = b'* cinv * b ;
  else
    b    = (X'*X)\(X'*y);
    k    = size(X,2);
    stat = ( b' * X' * X * b )./(p*(1-p));
  end
  dq_t = 1-chi2cdf(stat, k) ;
  catch
  dq_t = 0 ;
  end
end
 
