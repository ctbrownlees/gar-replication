function [status] = qr_oos( H , covs , createmfile )
% This function forecasts marginal and joint GaR using quantile regressions with
% several predictors.
%
% USAGE:
%   block_bootstrap(DATA,B,W)
%
% INPUTS:
%   H      - Vector of forecast horizons
%   covs   - Vector of coverage desired
%
% OUTPUTS:
%   BSDATA  - T by B matrix of bootstrapped data
%   INDICES - T by B matrix of locations of the original BSDATA=DATA(indexes);
%
% COMMENTS:
%   To generate bootstrap sequences for other uses, such as bootstrapping vector processes,
%   set DATA to (1:N)'.
%
% See also stationary_bootstrap

% Author: Andre B.M Souza and Christian T. Brownlees

%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO
%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clearing Workspace
%clear all
%clc
% For LASSO estimation
% for details, see
% https://faculty.fuqua.duke.edu/~abn5/belloni-software.html
%% Load Data Mfile, predictors and set params
load('../../../data/mfiles/data.mat')
covs  = 1 - covs  ;    % 1 - Coverage
% Load Predictors
sv = readtable('../../../data/clean/sv.csv','TreatAsEmpty',{'.','NA'});
ts = readtable('../../../data/clean/ts.csv','TreatAsEmpty',{'.','NA'});
hp = readtable('../../../data/clean/hp.csv','TreatAsEmpty',{'.','NA'});
cr = readtable('../../../data/clean/cr.csv','TreatAsEmpty',{'.','NA'});
% Collect dates
dates_sv = datenum(sv{:,2});
dates_ts = datenum(ts{:,2});
dates_hp = datenum(hp{sum(isnan(hp{:,3}))+1:end,2});
hp       = hp{sum(isnan(hp{:,3}))+1:end,3:end};
dates_cr = datenum(cr{:,2});
% Match with NFCI
sv_qr      = sv{ismember(dates_sv,dates_nfci),3};
ts_qr      = ts{ismember(dates_ts,dates_nfci),3:end};
hp_qr      = hp(ismember(dates_hp,dates_nfci),:);
cr_qr      = cr{ismember(dates_cr,dates_nfci),3:end};
% Compute the difference in sample size to keep track of indexes
dsmp = os_garch - os_nfci;
%% Start of the forecasting exercise
for h = H % For each horizon
  fprintf('Starting OOS QR for h=%s',h)
    for c = 1:length(covs) % for each coverage level
      rng(12)              % Set seed
      % First, clear all structs to make sure nothing remains from the last iteration
      clearing = whos();
      clearing1 = clearing(arrayfun(@(d) strcmp(d.class,'struct'),clearing));
      clearing2 = clearing(arrayfun(@(d) strcmp(d.class,'table'),clearing));
      clearing2(arrayfun(@(d) ismember(d.name,{'data','nfci_data'}),clearing2)) = [];

      if size(clearing2,1)>0
        clear (clearing2.name)
      end
      if size(clearing1,1)>0
        clear (clearing1.name)
      end
      % Set coverage level
      coverage = covs(c);
      % Out of Sample forecasting
      for t=os_garch:(T-h)
        [ ~, pc1QR , ~ ] = pca( gdp_garch( 1:t , : ) );
        % Use longest panel to compute factor
        factorQR = pc1QR( : , 1 );
        for j=1:N
          % Historical: Use whole sample
          Y  = gdp_garch( ( h + 1 ) : t , j );
          X  = ones( size( Y ) );
          b1 = rq( X , Y , coverage );
          M.QRUNC( t - os_garch + 1 , j ) =  b1 ;
          % Bonferroni
          b1 = rq( X , Y , coverage/N );
          M.QRBUNC(t - os_garch + 1,j) =  b1 ;

          % For all other models, use data only from 1973
          Y  = gdp_qr( ( h + 1 ):( t - dsmp ) , j );
          % Note that gdp_qr(t-dsmp,:) == gdp_garch(t , :)

          % Model 1: NFCI
          X  = [ones( size( Y ) ) , gdp_qr( 1:( t - h - dsmp ) , j ), ...
            nfci( 1:( t - h - dsmp ), j ) ] ;
          b1 = rq( X , Y , coverage);
          M.QRFIN( t - os_garch + 1 , j ) = [ 1 , gdp_qr( t - dsmp , j ) ,...
            nfci( t - dsmp , j ) ] * b1 ;
          % Bonferroni
          b1 = rq( X , Y ,coverage/N );
          M.QRBFIN( t - os_garch + 1 , j ) = [ 1 , gdp_qr( t - dsmp , j ) ,...
             nfci( t - dsmp , j ) ] * b1  ;

          % Model 2 NFCI TS
          X  = [ones( size( Y ) ) , gdp_qr( 1:(t - h - dsmp) , j ) , ...
            ts_qr( 1:t - h - dsmp , j ), nfci( 1:( t - h - dsmp ) , j )] ;
          b1 = rq( X , Y ,coverage);
          M.QRTS(t+1-os_garch,j) = [ 1 , gdp_qr(t - dsmp , j ) ,...
            ts_qr( t -dsmp , j ) , nfci( t - dsmp , j )] * b1 ;
          % Bonferroni
          b1 = rq( X , Y ,coverage/N);
          M.QRBTS(t+1-os_garch,j) = [ 1 , gdp_qr(t - dsmp , j ) ,...
            ts_qr( t -dsmp , j ) , nfci( t - dsmp , j )] * b1 ;

          % Model 3:  NFCI TS GF
          X  = [ones( size( Y ) ) , gdp_qr( 1:(t - h - dsmp ) , j ) ,...
            nfci( 1:( t - h - dsmp),j) ts_qr( 1:(t - h -dsmp) , j ), ...
            factorQR( (dsmp + 1) : ( t - h ) ) ] ;
          b1 = rq( X , Y ,coverage);
          M.QRGF( t - os_garch + 1 , j ) = [1 , gdp_qr( t - dsmp , j ) , ...
            nfci( t - dsmp , j ), ts_qr( t - dsmp , j ) , factorQR( t ) ] * b1 ;
          % Bonferroni
          b1 = rq( X , Y ,coverage / N );
          M.QRBGF( t - os_garch + 1 , j )  = [1 , gdp_qr(t - dsmp , j ) , ...
            nfci( t - dsmp , j ), ts_qr( t - dsmp , j) , factorQR( t ) ] * b1 ;

          % Full Model
          X  = [ ones( size( Y ) ) , gdp_qr( 1:( t - h - dsmp ) , j ) , ...
             factorQR( (dsmp + 1):( t - h) ) , nfci( 1:( t - h - dsmp ) , j ) ,...
             hp_qr( 1:( t - h - dsmp ) , j) , sv_qr( 1:( t - h -dsmp) ),...
             ts_qr( 1:( t - h - dsmp) , j)] ;
          b1 = rq( X , Y ,coverage);
          M.QRFULL( t + 1 - os_garch , j ) = [1 , gdp_qr( t - dsmp , j) ,...
             factorQR( t ) , nfci( t - dsmp , j ) , hp_qr( t - dsmp , j ) , ...
             sv_qr( t - dsmp) , ts_qr( t - dsmp , j )] * b1 ;
         % Bonferroni
          b1 = rq( X , Y ,coverage / N );
          M.QRBFULL( t + 1 - os_garch , j ) = [1 , gdp_qr( t - dsmp , j) ,...
             factorQR( t ) , nfci( t - dsmp , j ) , hp_qr( t - dsmp , j ) , ...
             sv_qr( t - dsmp) , ts_qr( t - dsmp , j )] * b1 ;

          % Lasso Quantiles
          Y    = gdp_qr( (h+1):( t - dsmp ) , j);
          X    = [ones( size( Y ) ) , gdp_qr(1:( t - h - dsmp) , j) ,...
             factorQR( (dsmp + 1):( t - h ) ) , nfci( 1:(t - h - dsmp) , j) ,...
             hp( 1:( t - h - dsmp ) , j) , sv_qr( 1:( t - h - dsmp) ),...
             ts_qr( 1:( t - h - dsmp),j) , cr_qr( 1:( t - h - dsmp) , j)] ;
          b1   = L1QR( Y , X , coverage);
          Xpos = X.*( abs( b1 ) > 10^(-6))';
          b1   = rq( Xpos , Y , coverage);
          M.QRLASSO( t + 1 - os_garch , j) =  [1 , gdp_qr( t - dsmp , j ) , ...
              factorQR( t ) , nfci( t - dsmp , j ) , hp_qr( t - dsmp , j ) ,...
              sv_qr( t - dsmp ) , ts_qr( t - dsmp , j ) ,cr_qr(t - dsmp , j ) ] * b1 ;
          % Bonferroni
          b1   = L1QR( Y , X , coverage / N );
          Xpos = X.*( abs( b1 )>10^(-6))';
          b1   = rq( Xpos , Y ,coverage / N);
          M.QRBLASSO(t + 1 - os_garch , j) = [1 , gdp_qr( t - dsmp , j ) , ...
              factorQR( t ) , nfci( t - dsmp , j ) , hp_qr( t - dsmp , j ) ,...
              sv_qr( t - dsmp ) , ts_qr( t - dsmp , j ) ,cr_qr( t - dsmp , j )] * b1 ;

        end
        fprintf(".")
    end
    fprintf(' QR OOS %s steps ahead done!' , h)
    save( sprintf( '../../../data/mfiles/qr/coverage_%.2f_%.0f_ahead.mat',...
                  1-coverage , h ) , 'M' ) % Save MFile to path
    saved(h) = isfile( sprintf( '../../../data/mfiles/qr/coverage_%.2f_%.0f_ahead.mat',...
                  1-coverage , h ) );
  end
end

if all(saved == 1)
    status = 1;
    fprintf('OOS QR done!')
else
    status = 0;
    fprintf('Problem at %s',find(saved==0,1));
end
