function [SLL,ht,zt,htfcast]=tarch_composite_LL(param,data,gdp_shocks,p,o,q,h)

% Initializing
[T , N ] = size(data);
ht = zeros(size(data));
htfcast = zeros([1 N]);
m = max([p,o,q,h]);
ht(1:m,:) = repmat(var(data),[h 1]);
% Forcing 0 if symmetric model.
if o<1
    param = [param(1) 0 param(2)];
end

for i = 1:N
  % Compute Variances
  for t = (m+1):T
   % ht(t,i) = ht(1,i).*(1-param(1)-param(3)-0.5*param(2))+...
  %            param(1)*data(t-h,i).^2 + ... % Note the direct GARCH
  %            param(2)*(data(t-h,i)<0)*data(t-h,i)^2+... 
  %            param(3)*ht(t-1,i);
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
%r1 = param(1) + param(2)*0.5 + param(3) >=1;
%r2 = (param(3)<0);
%r4 = ~(param(1)>0);
%r5 = ~(param(1)+param(2)>0);
%if o < 1
%    r3 = param(2) > 0;
%else
%    r3 = 0;
%end

%if any([r1 ,r2 ,r3 ,r4 , r5])
%  SLL=1000000;
%end

