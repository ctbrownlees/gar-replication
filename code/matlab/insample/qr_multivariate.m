function[status] = qr_multivariate( H , coverage )
global HOME
%% Load Data Mfile and set params
load( sprintf('%s/data/mfiles/data.mat', HOME ) )
sv = readtable( sprintf('%s/data/clean/sv.csv', HOME ),'TreatAsEmpty',{'.','NA'});
ts = readtable( sprintf('%s/data/clean/ts.csv', HOME ),'TreatAsEmpty',{'.','NA'});
hp = readtable( sprintf('%s/data/clean/hp.csv', HOME ),'TreatAsEmpty',{'.','NA'});
% Set dates
sv_dates = datenum(sv{:,2});
ts_dates = datenum(ts{:,2});
hp_dates = datenum(hp{:,2});

sv  = sv{ismember(sv_dates,dates_nfci),3};
ts  = ts{ismember(ts_dates,dates_nfci),3:end};
hp  = hp{ismember(hp_dates,dates_nfci),3:end};

% Scale vars
scale_var = @(x) (x - mean(x))./sqrt(var(x));
sv  = scale_var(sv);
ts  = scale_var(ts);
hp  = scale_var(hp);
nfci  = scale_var(nfci);

T     = size(gdp_qr,1);
B     = 1000;                      % Bootstrap sims
covs  = 1 - coverage ; % Desired Coverage = 1 - covs
%  Do we want to save tables?
createtables = 1;
createmfiles = 1;
%% Forecasting Exercise

for h = H
    for c = 1:length(covs)
      clear M
      rng(12)
      coverage = covs(c);
        [~, pc1QR , ~ ] = pca(gdp_qr);
        factorQR = pc1QR(:,1);

        Z = 1:T-h;
        Zb = block_bootstrap(Z',B,4);
        %  QREGS
        for j=1:N
          Y  = gdp_qr((h+1):T,j);
          % (1) NFCI only
          X  = [ones(size(Y)) ,gdp_qr(1:(T-h),j), nfci(1:(T-h),j)] ;
          b1 = rq( X , Y ,coverage);
          C1(j,:) = b1;

          for b=1:B
            bcoef(:,b) = rq( X(Zb(:,b),:) , Y(Zb(:,b)) ,coverage);
          end
          tstat_M1(:,j) = abs(C1(j,:)./std(bcoef'))>norminv(0.99);
          clear bcoef

          M.NFCI(:,j) =  X * b1 ;

          % (2) Real Variables
          X  = [ones(size(Y)) ,gdp_qr(1:(T-h),j),ts(1:(T-h),j) ,...
              nfci(1:(T-h),j) ] ;
          b1 = rq( X , Y ,coverage);
          C2(j,:) = b1;
          M.TS(:,j) =  X * b1 ;
          for b=1:B
            bcoef(:,b) = rq( X(Zb(:,b),:) , Y(Zb(:,b)) ,coverage);
          end
          tstat_M2(:,j) = abs(C2(j,:)./std(bcoef'))>norminv(0.99);
          clear bcoef
          % (3) TS + GF
          X  = [ones(size(Y)) ,gdp_qr(1:(T-h),j),factorQR(1:T-h), ts(1:(T-h),j) ,...
              nfci(1:(T-h),j) ] ;
          b1 = rq( X , Y ,coverage);
          C3(j,:) = b1;
          M.GFTS(:,j) =  X * b1 ;
          for b=1:B
            bcoef(:,b) = rq( X(Zb(:,b),:) , Y(Zb(:,b)) ,coverage);
          end
          tstat_M3(:,j) = abs(C3(j,:)./std(bcoef'))>norminv(0.99);
          clear bcoef

          % (4) Full Model
          X  = [ones(size(Y)) ,gdp_qr(1:(T-h),j), nfci(1:(T-h),j),ts(1:(T-h),j) , sv(1:T-h),...
              factorQR(1:T-h), hp(1:(T-h),j) ] ;
          b1 = rq( X , Y ,coverage);
          C4(j,:) = b1;
          M.FULL(:,j) =  X * b1 ;
          for b=1:B
            bcoef(:,b) = rq( X(Zb(:,b),:) , Y(Zb(:,b)) ,coverage);
          end
          tstat_M4(:,j) = abs(C4(j,:)./std(bcoef'))>norminv(0.99);
          clear bcoef
        end
        smaller = @(x) gdp_qr(h+1:end,:) < x ;           % Compute Hits
        H  = structfun(smaller , M , 'UniformOutput' , false );
    for i =1:N
           dqdirNFCI(i)   = dq_stat_direct(coverage, H.NFCI(:,i),h);
           dqdirREAL(i)  = dq_stat_direct(coverage, H.TS(:,i),h);
           dqdirFIN(i)  = dq_stat_direct(coverage, H.GFTS(:,i),h);
           dqdirFULL(i) = dq_stat_direct(coverage, H.FULL(:,i),h);
    end
      tick_loss   = @(x) mean( (gdp_qr(h+1:end,:)-x).*(coverage - ((gdp_qr(h+1:end,:)-x)<0) ) );
      Tloss     = structfun(tick_loss,M,'UniformOutput',false);
      %% Collecting variable specific coefficients.
      CY   = [C1(:,2)            C2(:,2)         C3(:,2)           C4(:,2)];
      CNFCI= [C1(:,3)            C2(:,end)       C3(:,end)         C4(:,3)];
      CTSP = [ zeros([N 1])      C2(:,3)         C3(:,4)           C4(:,4)];
      CFAC = [ zeros([N 1])      zeros([N 1])    C3(:,3)           C4(:,end-1)];
      CRV  = [ zeros([N 1])     zeros([N 1])    zeros([N 1])       C4(:,5)];
      CHP  = [ zeros([N 1])     zeros([N 1])    zeros([N 1])       C4(:,end)];

      %% And significances
      SY     = [ mean(tstat_M1(2,:)) mean(tstat_M2(2,:))       mean(tstat_M3(2,:))   mean(tstat_M4(2,:))];
      SNFCI  = [ mean(tstat_M1(3,:)) mean(tstat_M2(end,:))     mean(tstat_M3(end,:)) mean(tstat_M4(3,:))];
      STSP   = [ missing             mean(tstat_M2(end-1,:))   mean(tstat_M3(4,:))   mean(tstat_M4(4,:))];
      SFAC   = [ missing                 missing               mean(tstat_M3(3,:))   mean(tstat_M4(end-1,:))];
      SRV    = [ missing                 missing                  missing            mean(tstat_M4(5,:))];
      SHP    = [ missing                 missing                  missing            mean(tstat_M4(end,:))];

      %% Organizing quantiles
      AR   = [quantile(CY,0.25); quantile(CY,0.5); quantile(CY,0.75)];
      NFCI = [quantile(CNFCI,0.25); quantile(CNFCI,0.5); quantile(CNFCI,0.75)];
      FAC  = [quantile(CFAC,0.25); quantile(CFAC,0.5); quantile(CFAC,0.75)];
      TSP  = [quantile(CTSP,0.25); quantile(CTSP,0.5); quantile(CTSP,0.75)];
      CRV  = [quantile(CRV,0.25); quantile(CRV,0.5); quantile(CRV,0.75)];
      CHP  = [quantile(CHP,0.25); quantile(CHP,0.5); quantile(CHP,0.75)];
      DQ   = [mean(dqdirNFCI>0.01) mean(dqdirREAL>0.01) mean(dqdirFIN>0.01) mean(dqdirFULL>0.01)];
%% TABLE

QTB.AR   = AR';
QTB.TAR  = [missing  SY(1) missing  ; missing  SY(2) missing  ;missing  SY(3) missing ; missing  SY(4) missing ];
QTB.NFCI = NFCI';
QTB.TNFCI  = [ missing SNFCI(1) missing ;missing SNFCI(2) missing ; missing  SNFCI(3) missing ; missing  SNFCI(4) missing ];
QTB.TSP  = TSP' ;
QTB.TTSP  = [missing  STSP(1) missing  ; missing  STSP(2) missing  ; missing  STSP(3) missing  ; missing  STSP(4) missing ];
QTB.FAC  = FAC' ;
QTB.TFAC  = [missing  SFAC(1) missing  ; missing  SFAC(2) missing  ; missing  SFAC(3) missing  ; missing  SFAC(4) missing ];
QTB.CRV  = CRV' ;
QTB.TRV  = [missing  SRV(1) missing  ; missing  SRV(2) missing  ; missing  SRV(3) missing  ; missing  SRV(4) missing ];
QTB.CHP  = CHP';
QTB.TCHP  = [missing  SHP(1) missing   ; missing  SHP(2) missing  ; missing  SHP(3) missing  ;  missing  SHP(4) missing ];

TB{h}    = rows2vars(splitvars(struct2table(QTB)));
TB{h}    = [TB{h} ; ...
    table({'TL'}, mean(Tloss.NFCI) , mean(Tloss.TS) ,mean(Tloss.GFTS),mean(Tloss.FULL),'VariableNames',TB{h}.Properties.VariableNames) ; ...
    table({'DQ'}, DQ(1)*100 , DQ(2)*100 ,DQ(3)*100,DQ(4)*100,'VariableNames',TB{h}.Properties.VariableNames) ; ...
  ];
TB{h}.Properties.VariableNames = {'Predictors',sprintf('M1%0.f',h),sprintf('M2%0.f',h),sprintf('M3%0.f',h),sprintf('M4%.0f',h)};


        end
end
TBout =[ TB{1} TB{2}(:,2:end) TB{3}(:,2:end) TB{4}(:,2:end)];
writetable(TBout,sprintf('%s/tables/raw/insample/qr_models.csv' , HOME ),...
    'WriteRowNames',true)
[status] = isfile(sprintf('%s/tables/raw/insample/qr_models.csv' , HOME ));
if status == 1
  fprintf('QR Screening done \n')
else
  fprintf('Problem with QR Screening Table \n')
end
