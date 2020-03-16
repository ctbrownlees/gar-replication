function distributions_oos( H )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script constructs out of sample predictive densities as in
% Adrian et al (2019).
%
% INPUTS:
%   H     - Vector of forecast horizons
% 
% OUTPUTS:
%  This script creates predictive densities in .mat format in the folder
%  HOME/output 
%
% Authors: Christian Brownlees and Andre B.M. Souza
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covs = [0.05 0.25 0.75 0.95];
% We use the same quantiles as Adrian et al.
global HOME
rng(12)
for h = H
  % We need 4 quantiles to fit each skew T, so we group them below.
  for c = 1:length( covs )
    coverage = 1 - covs( c );
    % Data
    load( sprintf( '%s/data/input/data.mat' , HOME ))
    % QR results
    load( sprintf( '%s/data/output/QR_%.2f_%.0f_ahead.mat' , HOME , coverage , h ) )
    % Coverage vs 1-coverage
    % Drop the multivariate results
    marginal = ~contains(fieldnames(M),"B");
    qrNames  = fieldnames(M);
    qrNames  = qrNames(marginal);
    M    = struct2cell(M);
    M    = M(marginal);
    
    for i = 1:sum(marginal)
      QR{i,1}(:,:,c)     = M{i,1};
    end        

    % GARCH results
    load(sprintf('%s/data/output/GARCH_%.2f_%.0f_ahead.mat', HOME , coverage , h ) )
    % Drop the multivariate results
    marginal = contains(fieldnames(M),"marg");
    garchNames = fieldnames(M);
    garchNames = garchNames(marginal);
    M    = struct2cell(M);
    M    = M(marginal);
    for i = 1:sum(marginal)
      GARCH{i,1}(:,:,c) = M{i,1};
    end
  end
  clear M
  TMP = {QR{:} GARCH{:}};
  % Pasting all models together
  models = { qrNames{:} garchNames{:} };
  for i=1:(size(models,2))
    M.(models{i}) = TMP{i};
  end
  clear TMP
  
  % GDP against which we compare our forecasts:
  gdp_os      = gdp_garch( ( os_garch + h ) : T , :  );
  
  % Interpolate quantiles for each model, country and date and compute 
  % PITS and Predictive Sccores. This may take a
  % while...
  OS = size(QR{1,:},1);
  QQ = [0.05 0.25 0.75 0.95];
  PredScore = zeros(size(models,2),N,OS);
  PIT       = zeros(size(models,2),N,OS);
  LC = zeros(size(models,2),N,OS);
  SC = zeros(size(models,2),N,OS);
  SH = zeros(size(models,2),N,OS);
  DF = zeros(size(models,2),N,OS);
  for m=1:size(models,2)
  qt = M.(models{m});
    for cty = 1:N
      for t = 1:OS
        tmp = squeeze(qt(t,cty,:))';
        if t==1 
            [LC(m,cty,t),SC(m,cty,t),SH(m,cty,t),DF(m,cty,t)] = ...
                QuantilesInterpolation(sort(tmp),sort(QQ));
        else % provide initial guesses
            [LC(m,cty,t),SC(m,cty,t),SH(m,cty,t),DF(m,cty,t)] = ...
                QuantilesInterpolation(sort(tmp),sort(QQ),LC(m,cty,t-1),...
                SC(m,cty,t-1),SH(m,cty,t-1));
        end
         PredScore(m,cty,t) = dskt(gdp_os(t,1), LC(m,cty,t),...
           SC(m,cty,t), SH(m,cty,t), DF(m,cty,t));
         PIT(m,cty,t)  = pskt(gdp_os(t,cty), LC(m,cty,t),...
           SC(m,cty,t), SH(m,cty,t), DF(m,cty,t));
      end
     sprintf('Model: %.2f , Country : %.2f',m/size(models,2),cty/N)
    end
  end
  save(sprintf('data/output/distributions_%.0f_ahead.mat',h),'PredScore','PIT','LC','SC','SH','DF')
clear PredScore PIT LC SC SH DF  
  
end

 
fprintf('\nCSV files created, check the tables folder \n')

