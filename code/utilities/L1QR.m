%
% This is an initial implementation of L1QR 
% It relies on SDPT3 which needs to be downloaded and installed separetely.
% If lambda is omitted from the call, the choice of penalty level is done 
% as suggested in Belloni and Chernozhukov(2011).
% Each column is normalized.

function [ beta_L1QR ] = L1QR ( Y, X, tau, lambda )

NumSim = 5000;
eps    = 1.0e-10;
[ n p ] = size(X); 

if (nargin == 3)
    lambda_p    = L1QR_simulate_lambda ( X, tau, NumSim, n );
else
    lambda_p = lambda;
end
NormXX = zeros(p,1);
for j = 1 : p
    NormXX(j) = norm(X(:,j)/sqrt(n));  
end

%%% Parameters of L1-QR solver
OPTIONS.vers = 0; OPTIONS.gam = 0; OPTIONS.predcorr = 1; OPTIONS.expon = 1; OPTIONS.gaptol = 1.0e-10;
OPTIONS.inftol = 1.0e-10; OPTIONS.steptol = 1.0e-7; OPTIONS.maxit = 100; OPTIONS.printlevel = 0;
OPTIONS.stoplevel = 1; OPTIONS.scale_data = 0; OPTIONS.spdensity = 0.4; OPTIONS.rmdepconstr = 0; OPTIONS.cachesize = 128;
OPTIONS.smallblkdim = 15; OPTIONS.parbarrier = []; OPTIONS.schrfun = []; OPTIONS.schurfun_par = [];
%%% This version call SDPT3
[hn, opt_val_d, total_time, beta_L1QR, info] = L1QR_compact_formulation ( X , Y, tau, [ 0 ; lambda_p*( NormXX ) ; 0 ], OPTIONS.printlevel, OPTIONS, 0 );
%%
end

function [ lambda_final ] = L1QR_simulate_lambda ( XX, tau, NumSim, n )

NumSim = max(NumSim, n);

[ Numrows, NumColumns ] = size( XX );

U = rand(Numrows, NumSim );
NormXX = zeros(NumColumns,1);
lambda = zeros(NumSim,1);

for j = 1 : NumColumns
    NormXX(j) = norm(XX(:,j));
end    

for k = 1 : NumSim     
    lambda(k) = max( abs((XX'*( (U(:,k)<tau) - tau ) )./NormXX) );        
end

lambda_final = quantile(lambda, 0.9  ); 

end

function [opt_val_p, opt_val_d, total_time, beta, info] = L1QR_compact_formulation ( X, Y, tau, lambda_p, PrintLevel, OPTIONS, Mode )

start_time = cputime;

[ nX pX ] = size(X);

[ nY p1 ] = size(Y) ;

if ( abs(nX-nY) > 0 )
    fpritnf('Mismatch on the dimensions of X and Y\n');
    opt_val_p = 0;
    opt_val_d = 0;
    total_time = 0;
    beta = 0;
    info = 0;
    return;
end

n = nX;
p = pX;

K.f = p + p;
K.l = p + p + n + n ;
if ( max(size(lambda_p)) ==  1)
    lambda = lambda_p * ones( p, 1 );
else
    lambda = lambda_p;
end
%          beta            t         
c = [ zeros( p, 1) ; lambda ; tau * ones(n,1) ; (1-tau)*ones(n,1) ;  zeros( p + p, 1) ];
b = [ Y ; zeros( p + p , 1 ) ];
%     beta                 t                  u             v     (t_pos_slac,t_neg_slack)
%                                                                                     
%              
A = [  X               zeros(n,p)        eye(n,n)      -eye(n,n)       zeros(n,p+p) ;
     -eye(p,p)          eye(p,p)         zeros(p,n)     zeros(p,n)   -eye(p,p)  zeros(p,p) ;     
      eye(p,p)          eye(p,p)         zeros(p,n)     zeros(p,n)    zeros(p,p) -eye(p,p)  ];


pars.fid = 0;
if ( PrintLevel > 0 )
    pars.fid = 1;
end
pars.maxiter = 150;
pars.eps = 1.0e-8;
pars.bigeps = 1.0e-8;
pars.alg = 2;
pars.theta = 0.25;
pars.stepdif = 2;
pars.w = [1 1 ];
pars.numtol = 5 *1.0e-7;
pars.bignumtol = 0.9;
pars.numlvl = 0;
pars.denf =10;
pars.denq = 0.75;
pars.free = 1;
pars.sdp = 1;
pars.chol.skip = 1;
pars.chol.canceltol = 1.0e-12;
pars.chol.abstol = 1.0e-20;
pars.chol.maxuden = 500;
pars.cg.maxiter = 49;
pars.cg.refine = 1;
pars.cg.stagtol = 5*1.0e-14;
pars.cg.restol = 5*1.0e-3;
pars.cg.qprec = 1;
pars.vplot = 0;
pars.stopat = -1;
pars.errors = 0;


if ( Mode )
    [x,y,info] = sedumi(A,b,c,K,pars);
    opt_val_p = c'*x;
    opt_val_d = b'*y;

    beta = x(1:p);
else
  
    blk{1,1} = 'u';
    blk{1,2} = p+p;
    Avec{1} = A(:,1:blk{1,2});
    blk{2,1} = 'l';
    blk{2,2} = p+p+n+n;
    Avec{2} = A(:,blk{1,2}+1:blk{1,2}+blk{2,2});
    CC{1} = c(1:blk{1,2});
    CC{2} = c(blk{1,2}+1:blk{1,2}+blk{2,2});
    
    
    [obj,X,y,Z,info,runhist] = sqlp(blk,Avec,CC,b,OPTIONS);

    opt_val_p = obj(1);
    opt_val_d = obj(2);

    beta = X{1}(1:p);
end



total_time = cputime - start_time;

if ( PrintLevel > 0 )
    %fprintf('L1 norm Distance (within tolerance): %f    Total time: %f\n', opt_val_p, total_time);
end

end
