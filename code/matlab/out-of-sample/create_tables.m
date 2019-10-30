function [status] = create_tables( HAHEAD , coverages )
global HOME
rng(12)
B = 5000;
%% Body
for h = HAHEAD   
  for c=1:length(coverages)
    coverage=1-coverages(c);
    % Data
    load(sprintf('%s/data/mfiles/data.mat', HOME ))
    % QR results
    load(sprintf('%s/data/mfiles/qr/coverage_%.2f_%.0f_ahead.mat', HOME , 1 - coverage ,h))
    QR     = M;
    % GARCH results
    load(sprintf('%s/data/mfiles/garch/coverage_%.2f_%.0f_ahead.mat', HOME , 1 - coverage,h))
    % Joining into one struct
    fQ  = fieldnames(QR);     
    for m = 1: length(fQ)
      M.(fQ{m}) = QR.(fQ{m});
    end
    % Define out of sample variables:
    % GDP:
    gdp_os      = gdp_garch( ( os_garch + h ):T,: );
    % GDP agaisnt which we compare our forecasts;
    gdp_dq      = gdp_garch( os_garch:(T - h) , : );
    % GDP available when the forecast was made - for DQ.
    dates_os    = dates_garch( ( os_garch + h ): T , : );
    nfci_os     = nfci( os_nfci:end-h , : );
    % NFCI for DQ
    global_nfci = mean( nfci_os , 2 );
    %      [~, global_nfci , ~ ] = pca(nfci_os);
    %      global_nfci = global_nfci2(:,1);

    [~, global_real , ~ ] = pca( gdp_garch( os_garch:end-h , : ) );
    global_real = global_real( : , 1 );

    % Some helpful functions
    smaller      = @( x ) gdp_os <= x ;           % Compute Hits
    uni_hit      = @( x ) any( gdp_os <= x ,2 ) ; % Compute Uniform Hits
    mar_cov      = @( x ) mean( mean( x , 1 )'    ) ; % Marginal Coverage      
    len          = @( x ) sum( max( 0 , quantile( gdp_garch , 0.99 ) - x )) ...
    ./ (sum( quantile( gdp_garch , 0.99 ) >= x ) );
    dq_stat_cov  = @( x ) dq_stat_direct_hac_c( coverage , x , h );
    dq_fin_uni   = @( x ) aug_dq_stat_hac_c( coverage,x,global_nfci , h );
    dq_real_uni  = @( x ) aug_dq_stat_hac_c( coverage,x,global_real,h);
    uc_cov       = @( x ) cons_dq_stat_hac_c( coverage,x,h);
    tick_loss    = @( x ) mean( ( gdp_os - x ).*( coverage - ( ( gdp_os - x ) < 0 ) ) );
    % Hits        
    H  = structfun(smaller , M , 'UniformOutput' , false ); 
    % Uniform Hits
    UH = structfun(uni_hit , M , 'UniformOutput' , false );

    % Uniform Coverage
    C = table( 100 * ( 1 - structfun(@mean,UH)) , 'RowNames', fieldnames(H) , ...
    'VariableNames' , { 'UniformCoverage' } ) ;       

    % Marginal Coverage
    MC = table( 1 - structfun(mar_cov,H) , 'RowNames' , fieldnames(H), ...
    'VariableNames' , { 'MarginalCoverage' } );  

    % Length
    L    = structfun( len , M , 'UniformOutput' , false );
    avgL = table( structfun( @mean , L  ), 'RowNames' , fieldnames(H) , 'VariableNames' , {'Length'} ) ;
    L    = reshape(struct2array(L) , [] , numel( fieldnames( L ) ) ).' ;

    %% Table 4: Marginal Results
    % Collect marginal prediction regions
    margind = ~( ismember( fieldnames( M ),...
    { 'hjpr','bjprG','bjprT','bjprF','bonfG','bonfT','bonfF',...
    'QRBGF','QRBTS','QRBFULL','QRBFIN','QRBLASSO','QRBUNC' } ) ) ;
    fnames = fieldnames( M );
    marg   = fnames( margind );

    for k = 1:length(marg)
        for i=1:N
          % DQ UNC
          UUC(k,i)  = cons_dq_stat_hac_c(coverage, H.( marg{k} )( : , i ) , h );
          % Each country is backtested agaisnt its own data
          UDQ(k , i)  = dq_stat_direct_hac_c(coverage, H.( marg{k} )( : , i ) , h );
          UADQ(k , i) = aug_dq_stat_hac_c(coverage, H.( marg{k} )( : , i ), nfci_os( : , i ) , h );
          RADQ(k , i) = aug_dq_stat_hac_c(coverage, H.( marg{k} )( : , i ), gdp_dq( : , i ) , h );

          % Equal predictive ability tests: Unconditional
          lossh = ( gdp_os(: , i ) - M.QRUNC( : , i ) ).*( coverage - ( gdp_os( : , i ) - M.QRUNC( : , i ) < 0 ) );
          loss = ( gdp_os( : , i ) - M.( marg{k} )( :, i ) ).*( coverage - ( gdp_os( : , i )-M.( marg{k} )( : , i ) < 0 ) );
          lossdiff = loss - lossh ;
          int      = ones( length( lossdiff ) , 1  );
          try
            [~,~,~,covr] = olsnw( lossdiff , int , 0);
          catch
            covr = 1;
          end
          DMstath(k,i) = ( mean( lossdiff ) /sqrt( covr ) );
          DMh(k,i) = 1-normcdf( DMstath(k,i) );
          % Equal predictive ability tests Financial
          lossfin = ( gdp_os(:,i) - M.QRFIN( : , i ) ).*( coverage - ( gdp_os( : ,i)-M.QRFIN( : , i ) < 0) );
          lossdiff = loss - lossfin ;
          int      = ones(length(lossdiff),1);
          try
            [~,~,~,covr] = olsnw(lossdiff,int,0);
          catch
            covr = 1;
          end
          DMstatfin(k,i) = (mean(lossdiff)/sqrt(covr));
          DMfin(k,i) =  1-normcdf( DMstatfin(k,i) );
          % Equal predictive ability tests FGARCH
          lossf = ( gdp_os(:,i)-M.margF(:,i) ).*( coverage - (gdp_os(:,i)-M.margF(:,i)<0) );
          lossdiff = loss - lossf ;
          int      = ones(length(lossdiff),1);
          try
            [~,~,~,covr] = olsnw(lossdiff,int,0);
          catch
            covr = 1;
          end
          DMstatf(k,i) = (mean(lossdiff)/sqrt(covr));
          DMf(k,i) = 1-normcdf( (mean(lossdiff)/sqrt(covr)));
          % Equal predictive ability tests GARCH
          lossg = ( gdp_os(:,i)-M.margG(:,i) ).*( coverage - (gdp_os(:,i)-M.margG(:,i)<0) );
          lossdiff = loss - lossg ;
          int      = ones(length(lossdiff),1);
          try
            [~,~,~,covr] = olsnw(lossdiff,int,0);
          catch
            covr = 1;
          end
          DMstatg(k,i) = (mean(lossdiff)/sqrt(covr));
          DMg(k,i) =  1-normcdf( DMstatg(k,i) );
        end
      end

      U_MC     = table2array([ MC(margind,1)]);
      U_L      = table2array([ avgL(margind,1)]);

      Tloss     = structfun(tick_loss,M,'UniformOutput',false);
      Tloss     = reshape(struct2array(Tloss),[],numel(fieldnames(M)))';
      Tloss     = Tloss(margind,:);

      TB2.Mean    = [ 100*U_MC ];
      TB2.Length  = [U_L ];
      TB2.DQUNC    = [100*mean(UUC>0.01,2)];
      TB2.DQHITS    = [100*mean(UDQ>0.01,2)];
      TB2.DQFIN = [100*mean(UADQ>0.01,2)];
      TB2.DQREAL= [100*mean(RADQ>0.01,2)];
      TB2.Tloss   = [mean(Tloss,2)];

      TB2out  = struct2table(TB2);
      TB2out.Properties.RowNames = ["GARCH" , "GJR-GARCH" , "F-GARCH","Historical", "QR-NFCI" , "QR-NFCI+TS" ,...
      "QR-NFCI+TS+GF","Full","Lasso"]' ;
      TB2out = [TB2out(4,:) ; TB2out(5:end,:) ;TB2out(1:3,:)];

      writetable(TB2out,...
      sprintf('%s/tables/raw/Marginal_Coverage_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h),...
      'WriteRowNames',true)
      status(1) = isfile( sprintf('%s/tables/raw/Marginal_Coverage_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h));
      %% Table 5: Selected IMF Countries
      % IMF Countries
      % Countries for which IMF describes results
      imf = {'KOR' , 'AUS' , 'BRA' , 'CAN' , 'CHL' , 'CHI' , 'FRA' , 'DEU' ,...
      'IND' , 'IDN' , 'ITA' , 'JPN' , 'MEX' , 'RUS' , 'SAF' , 'ESP',...
      'SWE' , 'SWI' , 'TUR' , 'USA' , 'GBR' };

      names = data.Properties.VariableNames( 3:end );
      imf_names = names( ismember( names , imf ));
      imf_sel   = find( ismember( names , imf) );


      select      = ismember(fieldnames(M),{'QRUNC','QRFIN','margG'}) ;
      selectU     = ismember(fieldnames(M),{'QRUNC'}) ;
      selectFIN   = ismember(fieldnames(M),{'QRFIN'}) ;
      selectGARCH   = ismember(fieldnames(M),{'margG'}) ;
      select_marg = ismember(marg,{'QRUNC','QRFIN','margG'});
      select_unc  = ismember(marg,{'QRUNC'});
      select_fin  = ismember(marg,{'QRFIN'});
      select_gar  = ismember(marg,{'margG'});

      covIMF    = 100* (1- [ mean(H.QRUNC(:,imf_sel)) ; mean(H.QRFIN(:,imf_sel)) ; mean(H.margG(:,imf_sel)) ])';
      lenIMF    =  [ L(selectU,imf_sel) ; L(selectFIN,imf_sel) ; L(selectGARCH,imf_sel)]';
      TlossIMF  = [Tloss(select_unc,imf_sel);Tloss(select_fin,imf_sel);Tloss(select_gar,imf_sel)]';
      dmIMFh    = [ zeros(size(DMh(select_fin,imf_sel))) ; DMh(select_fin,imf_sel);DMh(select_gar,imf_sel)]';
      ucIMF     = [UUC(select_unc,imf_sel);UUC(select_fin,imf_sel);UUC(select_gar,imf_sel)]';
      dqIMF     = [UDQ(select_unc,imf_sel);UDQ(select_fin,imf_sel);UDQ(select_gar,imf_sel)]';
      dqfinIMF  = [UADQ(select_unc,imf_sel);UADQ(select_fin,imf_sel);UADQ(select_gar,imf_sel)]';
      dqrealIMF = [RADQ(select_unc,imf_sel);RADQ(select_fin,imf_sel);RADQ(select_gar,imf_sel)]';

      Historical    = [covIMF(:,1) lenIMF(:,1)  ucIMF(:,1) TlossIMF(:,1) dmIMFh(:,1)>0.99 dmIMFh(:,1)>0.95 dmIMFh(:,1)>0.90];
      GARCH   = [covIMF(:,3) lenIMF(:,3) ucIMF(:,3) TlossIMF(:,3) dmIMFh(:,3)>0.99 dmIMFh(:,3)>0.95 dmIMFh(:,3)>0.90];
      QRFIN   = [covIMF(:,2) lenIMF(:,2) ucIMF(:,2) TlossIMF(:,2) dmIMFh(:,2)>0.99 dmIMFh(:,2)>0.95 dmIMFh(:,2)>0.90];
      TB4 = table(Historical,QRFIN,GARCH);
      %TB4 = table(covIMF,lenIMF,ucIMF,TlossIMF);
      TB4.Properties.RowNames = imf_names;
      TB4 = sortrows(TB4,'RowNames');


      writetable(TB4,...
      sprintf('%s/tables/raw/IMF_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h),...
      'WriteRowNames',true)

      status(2) = isfile(sprintf('%s/tables/raw/IMF_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h)) ;

      %% Table 6: Uniform coverage
      DQHITS     = table(structfun(dq_stat_cov , UH ),'RowNames',fieldnames(H),'VariableNames',{'DQ'});
      DQFIN  = table(structfun(dq_fin_uni  , UH ),'RowNames',fieldnames(H),'VariableNames',{'DQFin'});
      DQREAL = table(structfun(dq_real_uni  , UH ),'RowNames',fieldnames(H),'VariableNames',{'DQReal'});
      DQUNC     = table(structfun(uc_cov , UH ),'RowNames' ,fieldnames(H),'VariableNames',{'UC'});

      TB = [ C  avgL  DQUNC DQHITS DQFIN DQREAL ];        

      % Choose rows that contain Uniform Prediction Regions
      unif_rows = ismember(TB.Properties.RowNames,...
      {'hjpr','bjprG','bjprT','bjprF','bonfG','bonfT','bonfF',...
      'QRBGF','QRBTS','QRBFULL','QRBFIN','QRBLASSO','margG','margT','margF'});
      TBout = TB(unif_rows,:);
      TBout.Properties.RowNames= ["Historical","GARCH + BJPR ","GARCH + Marg","GARCH + Bonf.","GJR-GARCH + BJPR", ...
      "GJR-GARCH + Marg","GJR-GARCH + Bonf.","F-GARCH + BJPR","F-GARCH + Marg","F-GARCH + Bonf.",...
      "QR-NFCI","QR-NFCI+TS","QR-NFCI+GF+TS","Full","Lasso"];

      % Rearranging rows
      TBout = [ TBout(1,:) ;TBout(11:end,:); TBout(2:3:9,:) ; TBout(4:3:11,:);...
      TBout(3:3:10,:)];

      % Writing
      writetable(TBout , ...
      sprintf('%s/tables/raw/Uniform_Coverage_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h),'WriteRowNames',true)
      status(3) = isfile( sprintf('%s/tables/raw/Uniform_Coverage_%s_%.0f_ahead.csv',...
      HOME , strrep(num2str(1-coverage,'%.2f'),'.',''),h));
      %% Table A2: Appendix
      TB3.TlossH     = [100*mean(DMh>0.99,2)];
      TB3.TlossFin   = [100*mean(DMfin>0.99,2)];
      TB3.TlossGARCH = [100*mean(DMg>0.99,2)];
      TB3.TlossFac   = [100*mean(DMf>0.99,2)];
      TB3out = struct2table(TB3);
      TB3out.Properties.RowNames = ["GARCH" , "GJR-GARCH" , "F-GARCH","Historical", "AR" , "Factor" ,...
      "Fin","Full","Lasso"]' ;
      TB3out = [TB3out(4,:) ;TB3out(5:end,:);TB3out(1:3,:)];

      writetable(TB3out,...
      sprintf('%s/tables/raw/Appendix_DM_%s_%.0f_ahead.csv',...
      HOME, strrep(num2str(1-coverage,'%.2f'),'.',''),h),...
      'WriteRowNames',true)
      status(4) = isfile(sprintf('%s/tables/raw/Appendix_DM_%s_%.0f_ahead.csv',...
      HOME, strrep(num2str(1-coverage,'%.2f'),'.',''),h));
    end
  end
