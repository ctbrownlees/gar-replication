function [status] = garch_insample( coverage )

global HOME

%% Load Data Mfile and set params
load(sprintf('%s/data/mfiles/data.mat' , HOME) )
B = 5000;                      % Bootstrap sims
coverage = 1 - coverage;       % Desired Coverage = 1 - covs

rng(12)
for i=1:N
  % Shocks to GDP growth
  X = [ones(T-1,1) gdp_garch(1:(T-1),i) ];
  Y = gdp_garch(2:T,i);
  b = (X'*X) \ (X'*Y) ;
  yhat(:,i) = X*b;
  phi_it(i,:) = b;
  % Forecast error
  gdp_shock(2:T,i)  = Y-X*b;
end

% Factor GARCH
[~, pc1 , ~ ]      = pca((gdp_shock(1:T,1:N)-nanmean(gdp_shock(1:T,1:N)))./sqrt(nanvar(gdp_shock(1:T,1:N))));
factor_shock       = pc1(:,1)/sqrt(nanvar(pc1(:,1)));
loadings_shock     = factor_shock(2:T)' * gdp_shock(2:T,1:N)/( factor_shock(2:T)'* factor_shock(2:T) );
idio_shock(:,1:N)  = gdp_shock(:,1:N) - factor_shock*loadings_shock;

[paramF , sigmaF, zFac , fcastF_t ] = tarch_composite(factor_shock(2:T),factor_shock(2:T),1,0,1,1);         

for i=1:N         
  % GARCH
  [paramG , sigmaGi, zGi ] = tarch_composite(gdp_shock(2:T,i),gdp_shock(2:T,i),1,0,1,1);
  pG(i,:) = paramG;
  sigmaG(:,i) = sigmaGi;
  zG(:,i)     = zGi;
  % TARCH
  [paramT , sigmaTi, zTi , fcastT_t ]  = tarch_composite(gdp_shock(2:T,i),gdp_shock(2:T,i),1,1,1,1); 
  pT(i,:)  = paramT;
  sigmaT(:,i) = sigmaTi;
  zT(:,i)     = zTi;
  % Factor GARCH
  [paramI , sigmaIi, zI , fcastI_t ]   = tarch_composite(idio_shock(2:T,i),idio_shock(2:T,i),1,0,1,1); 
  pI(i,:)=paramI;
  sigmaI(i,:) = sigmaIi;
end

%% Forecasting GaR

sigmaFG = ( loadings_shock.^2.*sigmaF + sigmaI');
zF      = gdp_shock(2:T,:)./sqrt(sigmaFG);

bdraws = datasample(1:(T-2),B);
bG   = zG(bdraws,:);
bT   = zT(bdraws,:);
bF   = zF(bdraws,:);

% GARCH 
M.margG  = yhat + sqrt( sigmaG )  .* quantile(bG,coverage,1);
% TARCH 
M.margT  = yhat + sqrt( sigmaT )  .* quantile(bT,coverage,1);
% Factor GARCH
M.margF  = yhat + sqrt( sigmaFG ) .* quantile(bF,coverage,1);

smaller = @(x) gdp_garch(2:end,:) < x ;           % Compute Hits
H  = structfun(smaller , M , 'UniformOutput' , false ); 

for i =1:N        
  archdY(i)        = arch_direct(gdp_shock(2:end,i),gdp_shock(2:end,i),1);
  archdG(i)        = arch_direct(zG(:,i),gdp_shock(2:end,i),1);
           archdT(i)        = arch_direct(zT(:,i),gdp_shock(2:end,i),1);
           archdF(i)        = arch_direct(zF(:,i),gdp_shock(2:end,i),1);
           dqdirG(i)        = dq_stat_direct(coverage, H.margG(:,i),1);
           dqdirT(i)        = dq_stat_direct(coverage, H.margT(:,i),1);
           dqdirF(i)        = dq_stat_direct(coverage, H.margF(:,i),1);
           
           lrG(i)           = 1- chi2cdf( 2*(tarch_composite_LL_zeros([ 0 0 ],gdp_shock(2:end,i),gdp_shock(2:end,i),1,0,1,1)-...
           tarch_composite_LL_zeros(pG(i,:),gdp_shock(2:end,i),gdp_shock(2:end,i),1,0,1,1)),1);
           lrT(i)           = 1- chi2cdf( 2*(tarch_composite_LL_zeros([ 0 0 0],gdp_shock(2:end,i),gdp_shock(2:end,i),1,1,1,1)-...
           tarch_composite_LL_zeros(pT(i,:),gdp_shock(2:end,i),gdp_shock(2:end,i),1,1,1,1)),2);
           lrF(i)           = 1- chi2cdf( 2*(tarch_composite_LL_zeros([ 0 0],idio_shock(2:end,i),idio_shock(2:end,i),1,0,1,1)-...
           tarch_composite_LL_zeros(pI(i,:),idio_shock(2:end,i),idio_shock(2:end,i),1,0,1,1)),1);  
           lrGammaT(i)      = 1- chi2cdf( 2*(tarch_composite_LL_zeros([ pG(i,1) 0 pG(i,2)],gdp_shock(2:end,i),gdp_shock(2:end,i),1,1,1,1)-...
           tarch_composite_LL_zeros(pT(i,:),gdp_shock(2:end,i),gdp_shock(2:end,i),1,1,1,1)),1);

end
    tick_loss = @(x) mean( (gdp_qr(2:end,:)-x(end-T_nfci+2:end,:)).*(coverage - ( (gdp_qr(2:end,:)-x(end-T_nfci+2:end,:) )<0) ) );
    Tloss     = structfun(tick_loss,M,'UniformOutput',false);

    %% GARCH VARS
    GARCH   = [quantile(pG,0.25); quantile(pG,0.5); quantile(pG,0.75)];
    GARCH_B = GARCH(:,2);
    GARCH_P = [quantile(pG(:,1)+pG(:,2),0.25); quantile(pG(:,1)+pG(:,2),0.5); quantile(pG(:,1)+pG(:,2),0.75)];
    TARCH   = [quantile(pT,0.25); quantile(pT,0.5); quantile(pT,0.75)];
    TARCH_B = TARCH(:,3);
    TARCH_G = TARCH(:,2);
    TARCH_P = [quantile(pT(:,1)+0.5*pT(:,2)+pT(:,3),0.25); quantile(pT(:,1)+0.5*pT(:,2)+pT(:,3),0.5);...
    quantile(pT(:,1)+0.5*pT(:,2)+pT(:,3),0.75)];
    FGARCH = [quantile(pI,0.25); quantile(pI,0.5); quantile(pI,0.75)];
    FGARCH_B = FGARCH(:,2);
    FGARCH_P = [quantile(pI(:,1)+pI(:,2),0.25); quantile(pI(:,1)+pI(:,2),0.5); quantile(pI(:,1)+pI(:,2),0.75)];
    PHI      = [quantile(phi_it(:,2),0.25) ; quantile(phi_it(:,2),0.5);quantile(phi_it(:,2),0.75)];
    C        = [quantile(mean(gdp_qr),0.25) ; quantile(mean(gdp_qr),0.5);quantile(mean(gdp_qr),0.75)];
    SIGMA    = [quantile(nanvar(gdp_shock),0.25);quantile(nanvar(gdp_shock),0.5);quantile(nanvar(gdp_shock),0.75)];
    ARCH     = [mean(archdG>0.01) mean(archdT>0.01) mean(archdF>0.01) ];
    TLG      = [quantile(Tloss.margG,0.25) quantile(Tloss.margG,0.5) quantile(Tloss.margG,0.75) ];
    TLT      = [quantile(Tloss.margT,0.25) quantile(Tloss.margT,0.5) quantile(Tloss.margT,0.75) ];
    TLFG     = [quantile(Tloss.margF,0.25) quantile(Tloss.margF,0.5) quantile(Tloss.margF,0.75) ];

    SKG      = [quantile(skewness(zG),0.25) quantile(skewness(zG),0.5) quantile(skewness(zG),0.75)];
    SKT      = [quantile(skewness(zT),0.25) quantile(skewness(zT),0.5) quantile(skewness(zF),0.75)];
    SKF      = [quantile(skewness(zF),0.25) quantile(skewness(zF),0.5) quantile(skewness(zT),0.75)];
    KURTG    = [quantile(kurtosis(zG),0.25) quantile(kurtosis(zG),0.5) quantile(kurtosis(zG),0.75)];
    KURTT    = [quantile(kurtosis(zT),0.25) quantile(kurtosis(zT),0.5) quantile(kurtosis(zF),0.75)];
    KURTF    = [quantile(kurtosis(zF),0.25) quantile(kurtosis(zF),0.5) quantile(kurtosis(zT),0.75)];
    DQ       = [mean(dqdirG>0.01) mean(dqdirT>0.01) mean(dqdirF>0.01)];
    LOAD     = [quantile(loadings_shock,0.25) quantile(loadings_shock,0.5) quantile(loadings_shock,0.75) ];
    %% GARCH TABLE
    GTB.Phi         = [ 0 0 0   PHI'  0 0 0 ]';
    GTB.Constant    = [ 0 0 0    C'   0 0 0 ]';
    GTB.Loadings    = [ 0 0 0  0 0 0   LOAD ]';
    GTB.Variance    = [ 0 0 0  SIGMA' 0 0 0 ]';

    GTB.Persistence = [GARCH_P' TARCH_P'   FGARCH_P']';
    GTB.Beta        = [GARCH_B' TARCH_B'   FGARCH_B']';
    GTB.Gamma       = [ 0 0 0   TARCH_G'     0 0 0  ]';
    GTB.Skewness    = [SKG         SKT        SKF   ]';
    GTB.Kurtosis    = [KURTG      KURTT      KURTF  ]';
    GTB.TickLoss    = [0 mean(Tloss.margG) 0 0 mean(Tloss.margT) 0 0 mean(Tloss.margF) 0]';
    GTB.LR          = 100*[0 1-mean(lrG>0.01) 0 0 1-mean(lrT>0.01) 0 0 1-mean(lrF>0.01) 0]';
    GTB.ARCH        = 100*[0 ARCH(1) 0 0 ARCH(2) 0 0 ARCH(3) 0]';
    GTB.DQ        = 100*[0 DQ(1) 0 0 DQ(2) 0 0 DQ(3) 0 ]';
    TB = struct2table(GTB);
    TB = rows2vars(TB);
    TB.Properties.VariableNames = {'Stats','GARCH025','GARCH050','GARCH075',...
    'TARCH025','TARCH050','TARCH075',...
    'FGARCH025','FGARCH050','FGARCH075'};
  
    writetable(TB,sprintf('%s/tables/raw/insample/garch_iterated_is.csv', HOME ), ...
        'WriteRowNames',true)
status = isfile(sprintf('%s/tables/raw/insample/garch_iterated_is.csv', HOME )) ;
if status == 1
  fprintf('GARCH in sample done! \n')
else
  fprintf('Problem in GARCH insample \n')
end
