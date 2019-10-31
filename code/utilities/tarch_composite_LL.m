function [SLL,ht,zt,htfcast]=tarch_composite_LL(param,data,gdp_shocks,p,o,q,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the composite likelihood of GARCH models, given 
% a set of parameters (param).
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
%   SLL         - Sum of Neg. Composite Likelihood
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
% Initializing
[T , N ] = size(data);
ht       = zeros(size(data));
htfcast  = zeros([1 N]);
m        = max([p,o,q,h]);
ht(1:m,:) = repmat(var(data),[h 1]);

% Forcing 0 if symmetric model.
if o<1
    param = [param(1) 0 param(2)];
end

for i = 1:N
  % Compute Variances
  for t = (m+1):T
    ht(t,i) = ht(1,i).*(1-param(1)-param(3)-0.5*param(2))+...
              param(1)*gdp_shocks(t-h,i).^2 + ... 
              param(2)*(gdp_shocks(t-h,i)<0).*gdp_shocks(t-h,i)^2+... 
              param(3)*ht(t-1,i);
   end
    % Compute LL
    [LL(i),~] = normloglik(data(:,i),0,ht(:,i)); 
    htfcast(i) = ht(1,i).*(1-param(1)-param(3)-0.5*param(2))+...
              param(1)*gdp_shocks(T+1-h,i).^2 + ... 
              param(2)*(gdp_shocks(T+1-h,i)<0).*gdp_shocks(T+1-h,i)^2+... 
              param(3)*ht(T,i);
 
end
zt = data./sqrt(ht);
SLL = sum(-LL);
% Restrictions on parameters (move this up)
r1 = param(1) + param(2)*0.5 + param(3) >=1;
r2 = (param(3)<0);
r4 = ~(param(1)>0);
r5 = ~(param(1)+param(2)>0);
if o < 1
    r3 = param(2) > 0;
else
    r3 = 0;
end

if any([r1 ,r2 ,r3 ,r4 , r5])
  SLL=1000000;
end


