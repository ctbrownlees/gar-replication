function [status] = garch_oos( H , covs )
global HOME
%% Load Data Mfile and set params
load(sprintf('data/mfiles/data.mat', HOME ))
B            = 5000         ;  % Number of Bootstrap sims
covs         = 1 - covs     ;  % Desired Coverage = 1 - covs
keep = whos();
keep = arrayfun(@(d) d.name,keep,'UniformOutput',false);
keep{end+1} = 'h';
keep{end+1} ='c';
keep{end+1} ='keep';
%% Forecasting Exercise
for h = H
  fprintf('Starting OOS GARCH for h=%s',h)
  for c = 1:length(covs)
    rng(12) % Set seed
    % First, clear everything from last iteration
    clearing = whos();
    clearing(arrayfun(@(d) ismember(d.name,keep),clearing)) = [];

    if size(clearing,1)>0
      clear (clearing.name)
    end

    coverage = covs(c);
    for t=is:(T-h)
      for i=1:N
        % Shocks to GDP growth and iterated forecast of cond.mean
        X = [ones(t-1,1) gdp_garch(1:(t-1),i) ];
        Y = gdp_garch(2:t,i);
        phi_it(i,:) = (X'*X) \ (X'*Y) ;
        gdp_shock(2:t,i)  = Y-X*phi_it(i,:)';
        b = [ phi_it(i,1)*( sum( phi_it(i,2).^((1:h)-1)) ) ; (phi_it(i,2)^h)];
        % Iterated Forecast of conditional mean
        yhat(t,i) = [1 gdp_garch(t,i)]*b;
      end

      %% Historical
      sgdp( t + 1 , : ) = ( gdp_garch(t,:) - mean(gdp_garch(1:t-1,:) )  ) ./ sqrt(var(gdp_garch(1:t-1,:),[],1) );

      % GARCH
      [paramG , sigmaG, zG, fcastG_t ] = tarch_composite(gdp_shock(h+1:t,:),gdp_shock(h+1:t,:),1,0,1,1);

      % TARCH
      [paramT , sigmaT, zT , fcastT_t ]= tarch_composite(gdp_shock(h+1:t,:),gdp_shock(h+1:t,:),1,1,1,1);

      % Factor GARCH (Decompose GDP shock into idio and common)
      beg = h+1; % Window size for loadings ( h+1 = recursive, t-x = X - window)

      [~, pc1 , ~ ] = pca((gdp_shock(1:t,1:N)-nanmean(gdp_shock(1:t,1:N)))./nanvar(gdp_shock(1:t,1:N)));
      factor = pc1(:,1);
      loadings(t,:)        = factor(beg:t)' * gdp_shock(beg:t,1:N)/( factor(beg:t)'* factor(beg:t) );
      idio_shock(1:t,1:N)  = gdp_shock(1:t,1:N) - factor(1:t)*loadings(t,:);
      idio_shock(1:t,N+1)  = factor;

      [paramF , sigmaF, zFac , fcastF_t ] = tarch_composite(factor(h+1:t),factor(h+1:t),1,0,1,1);
      [paramI , sigmaI, zI   , fcastI_t ] = tarch_composite(idio_shock(h+1:t,1:N),idio_shock(h+1:t,1:N),1,0,1,1);

      if t >= os_garch
        bdraws = datasample(1:(t-h),B);
        % Historical JPR
        zH    = sgdp(bdraws,:);
        dmax  = quantile( min(zH,[],2),coverage,1)';
        M.hjpr(t-os_garch+1,:) = mean(gdp_garch(1:t-1,:))' + sqrt(var( ( gdp_garch(1:t-1,:) ) ))'.*dmax;

        %% Iterated GARCH

        % Variances for variance targetting
        vtarget  = var(gdp_shock);
        vtargetI = var(idio_shock(:,1:N));
        vtargetF = var(factor);

        for q=1:h
          % Bootstrapping days
          bG   = zG(bdraws,:);
          bT   = zT(bdraws,:);
          bI   = zI(bdraws,:);
          bF   = zFac(bdraws,:);

          % Simulated paths
          % If h=1, use last gdp to forecast. Else, use previous
          % path.
          if q==1
            % Variance forecast from q-1 to q
            fcastG(q,:,:) = repmat(fcastG_t,[B 1]);
            fcastT(q,:,:) = repmat(fcastT_t,[B 1]);
            fcastI(q,:,:) = repmat(fcastI_t,[B 1]);
            fcastF(q,:,:) = repmat(fcastF_t,[B 1]);
            % h=1 Bootstrap returns
            ret_g(q,:,:)  = sqrt(fcastG_t).*bG;
            ret_t(q,:,:)  = sqrt(fcastT_t).*bT;
            ret_f(q,:,:)  = sqrt(fcastF_t).*bF;
            ret_i(q,:,:)  = sqrt(fcastI_t).*bI;
            ret_fg(q,:,:) = [loadings(t,:)'*ret_f(q,:,:)]' + squeeze(ret_i(q,:,:));
            % Path forecast
            path_g(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*gdp_garch(t,:) + squeeze(ret_g(q,:,:));
            path_t(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*gdp_garch(t,:) + squeeze(ret_t(q,:,:));
            path_f(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*gdp_garch(t,:) + squeeze(ret_fg(q,:,:));
          else
            % h>1 Bootstrap returns
            ret_g(q,:,:)  = squeeze(sqrt(fcastG(q,:,:))).*bG;
            ret_t(q,:,:)  = squeeze(sqrt(fcastT(q,:,:))).*bT;
            ret_f(q,:,:)  = sqrt(fcastF(q,:,:))'.*bF;
            ret_i(q,:,:)  = squeeze(sqrt(fcastI(q,:,:))).*bI;
            ret_fg(q,:,:) = [loadings(t,:)'.*ret_f(q,:,:)]' + squeeze(ret_i(q,:,:));

            % Path forecast
            path_g(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*squeeze(path_g(q-1,:,:)) + squeeze(ret_g(q,:,:));
            path_t(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*squeeze(path_t(q-1,:,:)) + squeeze(ret_t(q,:,:));
            path_f(q,:,:) = phi_it(:,1)' + phi_it(:,2)'.*squeeze(path_f(q-1,:,:)) + squeeze(ret_fg(q,:,:));
          end

          % q+1 variance forecast
          fcastG(q+1,:,:) = vtarget*(1-paramG(1) - paramG(2)) + ...
          paramG(1).*squeeze(ret_g(q,:,:)).^2 + paramG(2).*squeeze(fcastG(q,:,:));

          fcastT(q+1,:,:) = vtarget*(1- paramT(1) - paramT(2)*0.5 - paramT(3)) + paramT(1).*squeeze(ret_t(q,:,:)).^2 +...
          paramT(2).*squeeze(ret_t(q,:,:)).^2.*(squeeze(ret_t(q,:,:))<0) + paramT(3).*squeeze(fcastT(q,:,:));

          fcastI(q+1,:,:)   = vtargetI*(1-paramI(1)-paramI(2)) + ...
          paramI(1).*squeeze(ret_i(q,:,:)).^2 + paramI(2).*squeeze(fcastI(q,:,:));

          fcastF(q+1,:,:)   = vtargetF*(1-paramF(1)-paramF(2)) + ...
          paramF(1).*squeeze(ret_f(q,:,:)).^2 + paramF(2).*squeeze(fcastF(q,:,:));

          % Bootstrap new days
          bdraws = datasample(1:t-h,B);
        end

        % GARCH
        spaths = squeeze( (path_g(end,:,:) - mean(path_g(end,:,:)))./sqrt(var(path_g(end,:,:))));
        dmaxG   = quantile( min(spaths,[],2), coverage,1)';
        M.bjprG(t-os_garch+1,:) = mean(squeeze(path_g(end,:,:)))' + sqrt(var(squeeze(path_g(end,:,:))))'.*dmaxG;
        M.margG(t-os_garch+1,:) = quantile( squeeze( path_g(end,:,:)) , coverage, 1);
        M.bonfG(t-os_garch+1,:) = quantile( squeeze( path_g(end,:,:)) , coverage/N, 1);

        % TARCH
        spaths = squeeze( (path_t(end,:,:) - mean(path_t(end,:,:)))./sqrt(var(path_t(end,:,:))) );
        dmaxT   = quantile( min(spaths,[],2),coverage,1)';
        M.bjprT(t-os_garch+1,:) = mean(squeeze(path_t(end,:,:)))' + sqrt(var(squeeze(path_t(end,:,:))))'.*dmaxT;
        M.margT(t-os_garch+1,:) = quantile(  squeeze(path_t(end,:,:)) , coverage, 1);
        M.bonfT(t-os_garch+1,:) = quantile(  squeeze(path_t(end,:,:)) , coverage/N, 1);

        % Factor GARCH
        spaths = squeeze( (path_f(end,:,:) - mean(path_f(end,:,:)))./sqrt(var(path_f(end,:,:))));
        dmaxF   = quantile( min(spaths,[],2),coverage,1)';
        M.bjprF(t-os_garch+1,:) = mean(squeeze(path_f(end,:,:)))' + sqrt(var(squeeze(path_f(end,:,:))))'.*dmaxF;
        M.margF(t-os_garch+1,:) = quantile( squeeze( path_f(end,:,:)) , coverage, 1);
        M.bonfF(t-os_garch+1,:) = quantile( squeeze( path_f(end,:,:)) , coverage/N, 1);

      end
      fprintf('.')
    end
    fprintf(' GARCH OOS %s steps ahead done!' , h) 
    save(sprintf('%s/data/mfiles/garch/coverage_%.2f_%.0f_ahead.mat',HOME , ...
    1 - coverage , h ) , 'M' )
    status(h,c) = isfile(sprintf('%s/data/mfiles/garch/coverage_%.2f_%.0f_ahead.mat',HOME , ...
    1 - coverage , h ) );
  end
end
