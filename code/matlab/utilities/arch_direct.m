function[arch_t] = arch_direct(y,x,nahead)
  n    = length(x);
  y    = y(4+nahead:n).^2;
  X    = double([ones(size(y)), x(4:(n-nahead)).^2 , x(3:(n-1-nahead)).^2 , x(2:(n-2-nahead)).^2 , x(1:(n-3-nahead)).^2]);
  b    = (X'*X)\(X'*y);
  k    =  size(X,2);
  r2   = ((X*b-mean(X*b))'*(X*b-mean(X*b)))./((y-mean(y))'*(y-mean(y)));
  arch_t = 1-chi2cdf((n-4).*r2, k-1) ;
end
 