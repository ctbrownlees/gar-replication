function panel_qr_oos( H , covs )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function forecasts marginal and joint GaR using quantile regressions
% with several predictors.
%
% INPUTS:
%   H      - Vector of forecast horizons
%   covs   - Vector of coverage levels
%
% OUTPUTS:
%  This function creates .mat files in data/output. Each file contains
%  a struct, M, that contains the collection of forecasts for all models,
%  given a coverage level and a forecast horizon.
%
% SEE ALSO:
%  L1QR.m , rq.m
% COMMENTS:
%
% Authors: Christian Brownlees and Andre B.M. Souza
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global HOME
%% Load Data Mfile, predictors and set params
load( sprintf( '%s/data/input/data.mat' , HOME ))
covs  = 1 - covs  ;    % 1 - Coverage
% Load Predictors
sv = readtable( sprintf('%s/data/input/sv.csv' , HOME ),'TreatAsEmpty',{'.','NA'});
ts = readtable( sprintf('%s/data/input/ts.csv' , HOME ),'TreatAsEmpty',{'.','NA'});
hp = readtable( sprintf('%s/data/input/hp.csv' , HOME ),'TreatAsEmpty',{'.','NA'});
cr = readtable( sprintf('%s/data/input/cr.csv' , HOME ),'TreatAsEmpty',{'.','NA'});
% Collect dates
dates_sv = datenum( sv{ : , 2 });
dates_ts = datenum( ts{ : , 2 });
dates_hp = datenum( hp{ sum( isnan( hp{ : , 3 }))+1 : end , 2 });
hp       = hp{ sum( isnan( hp{ : , 3 }))+1 : end , 3 : end};
dates_cr = datenum( cr{ : , 2 });
% Match with NFCI
sv_qr      = repmat( sv{ ismember( dates_sv , dates_nfci), 3 }, 1,24);
ts_qr      = ( ts{ ismember( dates_ts , dates_nfci), 3 : end} );
hp_qr      = ( hp( ismember( dates_hp , dates_nfci), : )) ;
cr_qr      = ( cr{ ismember( dates_cr , dates_nfci) , 3 : end});
% Compute the difference in sample size to keep track of indexes
dsmp = os_garch - os_nfci;
keep = whos();
keep = arrayfun(@(d) d.name , keep , 'UniformOutput' , false);
keep{end+1} = 'h';
keep{end+1} ='c';
keep{end+1} ='keep';
% keep all these variables after each loop iteration, but remove
% everything else. 
%%%%%%%%%%%%%%%%%%%%% Forecasting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for h = H 
  fprintf('Starting OOS QR for h = %d\n',h)
  for c = 1:length(covs) % for each coverage level
    rng(12)              % Set seed
    % First, clear everything from last iteration
    clearing = whos();
    clearing( arrayfun(@(d) ismember(d.name,keep) , clearing)) = [];

    if size(clearing,1)>0
      clear (clearing.name)
    end

    % Set coverage level
    coverage = covs(c);
    % Out of Sample forecasting
    b1 = [];
    b2 = [];
    b3 = [];
    b4 = [];
    for t=os_garch:(T-h)
      [ ~, pc1QR , ~ ] = pca( gdp_garch( 1:t , : ) );
      % Use longest panel to compute factor
      factorQR = pc1QR( : , 1 );
      factorQR = repmat(factorQR,1,24);
      % For all other models, use data only from 1973
      Y  = gdp_qr( ( h + 1 ):( t - dsmp ) , : );
      % Note that gdp_qr(t-dsmp,:) == gdp_garch(t , :)

      % Model 1: NFCI
      X  = [ gdp_qr( 1:( t - h - dsmp ), :  ), nfci( 1:( t - h - dsmp ), : ) ] ;
      b1 = panel_rq( X , Y , coverage,b1);
      M.QRFIN( t - os_garch + 1 , : ) = [ b1(1:N) + ...
        + b1( (N+1): (2*N)).* gdp_qr( t - dsmp , : ) + nfci( t - dsmp,: ) * b1(end) ] ;
        
      % Bonferroni
      % b1 = rq( X , Y ,coverage/N );
      %M.QRBFIN( t - os_garch + 1 , : ) = [ b1(1:N) + ...
      %  + b1( (N+1): (2*N)).* gdp_qr( t - dsmp , : ) + nfci( t - dsmp ) * b1(end) ]  ;

      % Model 2 NFCI TS
      X  = [ gdp_qr( 1:(t - h - dsmp) , : ) , ts_qr( 1:t - h - dsmp ,:),...
        nfci( 1:( t - h - dsmp ), : )] ;
      b2 = panel_rq( X , Y ,coverage,b2);
      M.QRTS(t + 1 -os_garch , :) = [ b2(1:N) + ...
        + b2( (N+1): (2*N)).* gdp_qr( t - dsmp , : ) + ...
        + ts_qr( t - dsmp,: ) *b2(end-1) + nfci( t - dsmp,: ) * b2(end) ];
      % Bonferroni
      %b1 = rq( X , Y ,coverage/N);
      %M.QRBTS(t + 1 - os_garch , j) = [ 1 , gdp_qr(t - dsmp , j ) ,...
      %ts_qr( t - dsmp , j ) , nfci( t - dsmp , j )] * b1 ;
      % Model 3:  NFCI TS GF
      X  = [gdp_qr( 1:(t - h - dsmp ) , : ) ,...
      ts_qr( 1:(t - h - dsmp), : ), nfci( 1:( t - h - dsmp),:), ...
      factorQR( (dsmp + 1) : ( t - h ) , :  ) ] ;
      b3 = panel_rq( X , Y ,coverage,[]);
      M.QRGF( t - os_garch + 1 , : ) = [ b3(1:N) + ...
        + b3( (N+1): (2*N)).* gdp_qr( t - dsmp , : ) + ...
        + ts_qr( t - dsmp ) *b3(end-2) + nfci( t - dsmp ) * b3(end-1) + ...
        + factorQR(t,:) * b3(end)];
      % Bonferroni
      %b1 = rq( X , Y ,coverage / N );
      %M.QRBGF( t - os_garch + 1 , j )  = [1 , gdp_qr(t - dsmp , j ) , ...
      %nfci( t - dsmp , j ), ts_qr( t - dsmp , j) , factorQR( t ) ] * b1 ;
       % Full Model
      X  = [ gdp_qr( 1:( t - h - dsmp ) , : ) , ...
      ts_qr( 1:( t - h - dsmp),: ), nfci( 1:( t - h - dsmp ),: ),...
      factorQR( (dsmp + 1):( t - h),: ) , hp_qr( 1:( t - h - dsmp ),: ) , ...
      sv_qr( 1:( t - h -dsmp) ,: ) ] ;
      b4 = panel_rq( X , Y ,coverage,b4);
      M.QRFULL( t + 1 - os_garch , : ) = [ b4(1:N) + ...
        + b4( (N+1): (2*N)).* gdp_qr( t - dsmp , : ) + ...
        + ts_qr( t - dsmp ,: ) *b4(end-4) + nfci( t - dsmp , : ) * b4(end-3) + ...
        + factorQR(t, : ) * b4(end-2) + hp_qr( t - dsmp, : ) * b4(end-1) + ...
        + sv_qr( t - dsmp , :) * b4(end) ];
     
      % Bonferroni
      %b1 = rq( X , Y ,coverage / N );
      %M.QRBFULL( t + 1 - os_garch , j ) = [1 , gdp_qr( t - dsmp , j) ,...
      %factorQR( t ) , nfci( t - dsmp , j ) , hp_qr( t - dsmp , j ) , ...
      %sv_qr( t - dsmp) , ts_qr( t - dsmp , j )] * b1 ;
     fprintf('.')
  end
    fprintf('\n Panel - QR OOS %d-step-ahead done!\n' , h)
    save( sprintf( '%s/data/output/panel_QR_%.2f_%.0f_ahead.mat',...
    HOME , 1-coverage , h ) , 'M' ) % Save MFile to path
  end
end
