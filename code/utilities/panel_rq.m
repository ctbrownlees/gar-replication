function [p]=panel_rq(x,y,tau,guess);
% Panel Quantile Regression
% 
% USAGE: [p,stats]=quantreg(x,y,tau[,order,nboot]);
% 
% INPUTS: 
%   x,y: data that is fitted. (x and y should be columns)
%   First Column of X is each countrys lagged gdp, each country is allowed
%   to load individually on this column.
%   tau: quantile used in regression. 
%   order: polynomial order. (default=1)
%          (negative orders are interpreted as zero intercept)
%   nboot: number of bootstrap surrogates used in statistical inference.(default=200)
%

if isempty(guess)
  N = size(y,2);
  P = size(x,2)/N;
  for i=1:N
    tmp       = [x(:, (0:(P-1))*N + i)];
    pmean(i,:)= tmp\y(:,i);
  end
  guess = [ones(1,N) pmean(:,1)' mean(pmean(:,2:P))];
end

opt  = struct('MaxFunEvals',50000,'MaxIter',50000);
p=fminsearch(@(par)TickLoss(x,y,tau,par),guess,opt);
end

function [loss] = TickLoss(x,y,tau,theta)
N = size(y,2);
P = size(x,2)/N;
rho=@(r)sum(abs(r.*(tau-(r<0))));

for i=1:N
  L(i) = rho(y(:,i) - theta(i) - theta(N + i).*x(:,i) - ...
    x(:,((1:P-1)*N) + i)*theta((end-P+2):end)' ); % 1 is intercept, 2 is AR.
end
loss = sum(L);
end

