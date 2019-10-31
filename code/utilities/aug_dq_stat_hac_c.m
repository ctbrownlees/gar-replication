function[dq_t] = aug_dq_stat_hac_c(p,h,x,nahead)  
  n    = length(h);
  y    = h(5:n)-p;
  X    = double([ones(size(y)), x(4:(n-1)) , x(3:(n-2)) , x(2:(n-3)) , x(1:(n-4))]);
  k    = size(X,2);  
 try
   if nahead > 1
    [b , ~ ,~ , covmat] = olsnw_dq(y,X,0,nahead-1,p);
    cinv   = inv(covmat);
    stat = b'* cinv * b;
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
 
