function [status] = qr_univariate(x,xdates, H , covs , name)
global HOME

load(sprintf('%s/data/mfiles/data.mat',HOME))
B = 1000;
dates_keep = intersect(dates_nfci,xdates);
if isempty(x)
  dates_keep = dates_nfci;
end
gdp_keep   = gdp_qr(ismember(dates_nfci,dates_keep),:);
x          = x(ismember(xdates,dates_keep),:);

T = size(gdp_keep,1);

keep = whos();
keep = arrayfun(@(d) d.name,keep,'UniformOutput',false);
keep{end+1} = 'h';
keep{end+1} = 'TB';
keep{end+1} ='c';
keep{end+1} ='keep';
nfci = (nfci  - mean(nfci) ) ./ sqrt(var(nfci));


%% Benchmark

for h = H
    for c = 1:length(covs)

        rng(12)
        % First, clear everything from last iteration
        clearing = whos();
        clearing(arrayfun(@(d) ismember(d.name,keep),clearing)) = [];

        if size(clearing,1)>0
          clear (clearing.name)
      end
      coverage = 1 - covs(c);
      if name=='gf'
        [~, pc1QR , ~ ] = pca(gdp_qr);
        factorQR = pc1QR(:,1);
      end

      Z = 1:T-h;
      Zb = block_bootstrap(Z',B,4);
      %  QREGS
      for j=1:N
          Y  = gdp_keep((h+1):T,j);
          if name == 'gf'
            X  = [ones(size(Y))  , gdp_keep(1:(T-h),j) , factorQR(1:T-h)  ] ;
          elseif name == 'nfci'
            X  = [ones(size(Y))  , gdp_keep(1:(T-h),j) , nfci(1:T-h,j)  ] ;
          elseif size(x,2) > 1
            X  = [ones(size(Y))  , gdp_keep(1:(T-h),j) , x(1:T-h,j)  ] ;
          else
            X  = [ones(size(Y))  , gdp_keep(1:(T-h),j) , x(1:T-h)  ] ;
          end
          b1 = rq( X , Y ,coverage);
          COEF(j) = b1(end);
          M(:,j) = X * b1 ;
          clear bcoef
          for i=1:B
             bcoef(:,i) = rq( X(Zb(:,i),:) , Y(Zb(:,i)) ,coverage);
          end
          TSTAT(j) = b1(end)'./std(bcoef(end,:));

          Y  = gdp_keep((h+1):T,j);
          X  = [ones(size(Y))  , gdp_keep(1:(T-h),j) , nfci(1:T-h,j)  ] ;
          b1 = rq( X , Y ,coverage);
          BM(:,j) = X * b1 ;


        end
        TL  =  (mean( (gdp_keep(h+1:end,:)-M).*(coverage - ((gdp_keep(h+1:end,:)-M)<0) ) ));
        TL  = 1 - TL /(mean( (gdp_keep(h+1:end,:)-BM).*(coverage - ((gdp_keep(h+1:end,:)-BM)<0) ) ));

        coef = [quantile(COEF,0.25) quantile(COEF,0.5) quantile(COEF,0.75)];
        ts = [quantile(TSTAT,0.25) quantile(TSTAT,0.5) quantile(TSTAT,0.75)];

        TB{h} = rows2vars(table(coef',[0 mean(TL) 0]',100*[0 mean(abs(TSTAT)>norminv(0.99)) 0 ]'));
        TB{h}.Properties.RowNames = {'Coef','TickLoss','pass'};
        TB{h}.Properties.VariableNames = {'STATS',sprintf('Q025%.0f',h),sprintf('Q050%.0f',h),sprintf('Q075%.0f',h)};
      end
end

date_av = cellstr(repmat(sprintf('%.0fQ%.0f to %.0fQ%.0f',year(dates_keep(1)),quarter(dates_keep(1)),...
year(dates_keep(end)),quarter(dates_keep(end))),[3 1]));

TBout = [table(date_av) TB{1}(:,2:end) TB{2}(:,2:end) TB{3}(:,2:end) TB{4}(:,2:end)];
if exist('name','var')
     writetable(TBout,...
     sprintf('%s/tables/raw/insample/%s.csv', HOME , name ),'WriteRowNames',true)
end
status = isfile(sprintf('%s/tables/raw/insample/%s.csv', HOME ,name));
if status ==1
   fprintf('Univariate QR for %s done! \n',name)
 else
  fprintf('Problem in Univariate QR for %s \n',name)
  end
