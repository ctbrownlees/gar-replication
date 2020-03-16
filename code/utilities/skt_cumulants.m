function cum=skt_cumulants(location, scale, shape, df, n)
%skt_cumulants
%cumulants of the skew-t distribution
% 
%DESCRIPTION
% 
%cumulants of the skew-t distribution
% 
%USAGE
% 
%skt_cumulants(location, scale, shape, n)
% 
%ARGUMENTS
% 
%location   a vector of location parameters (default is 0).
%scale      a vector of scale parameters (default is 1).
%shape      a vector of shape parameter (default is 0).
%df         degrees of freedom (scalar); (default is df=Inf which corresponds to the
%           skew-normal distribution).
%n	        a scalar integer of the maximal order or cumulants required; it must be
%           from 1 to 4 and smaller than df (default is 4).
%cum        a vector of 4 elements which are taken to represent the first 4 cumulants
%           of a skew-t distribution (hence the second term must be positive)
%abstol     a scalar which regulates the accuracy of the cumulants matching (default 
%           value 1e-08).
% 
%VALUE
% 
%skt_cumulants computes the cumulants up to order n of the skew-t distribution with
%the selected parameters. The returned object is a vector of lemgth n if the parameters 
%are all scalar, otherwise a matrix with n coilumns. 
% 
%DETAILS
% 
%Expressions of the moments and other details on the skew-t distribution are given in 
%reference bellow. These formulae are used by skt_cumulants to compute the cumulants.
%
%skt_cumulants_inversion searches the set of shape and df parameters of the skew-t family,
%attempting to match the third and fourth cumulants with those of the supplied vector cum.
%This search is done numerically twice, once using optim and a second time using nlminb, to
%the accuracy abstol; the best matching solution is retained. If the required accuracy of the 
%matching is not achieved by any of the two methods, a warning message is issued. After this
%step, the other two parameters (location and scale) are computed via simple algebra.
%
%The moment generating function (hence the cumulant generating function)
%of the distribution is given in the refence below. 
% 
%REFERENCES
% 
%Azzalini, A. (1985). A class of distributions which includes the normal
%ones. Scand. J. Statist. 12, 171-178.
% 
%NOTE
%
%The joint use skt_cumulants_inversion and sample_centralmoments allows to fit a skew-t
%distribution by the method of moments; see example bellow. Note however, that for 
%stability reasons, this is not adopted as the standard method for producing initial values 
%of the MLE search.

%SEE ALSO
% 
%sn_cumulants, dskt, sample_centralmoments, optim, nlminb
% 
%EXAMPLES
% 
%a   = skt_cumulants(10, 2, -8, 5.2)    #Result is: 8.1310    3.0070  -12.2242  144.2662
%b   = skt_cumulants( 0, 1,  0, 5)      #Result is:     0    1.6667         0   16.6667
%c   = skt_cumulants( 0, 1,  3, 5)      #Result is: 0.9003    0.8561    1.6846   11.9140
%d   = skt_cumulants( 0, 1,  9, 5)      #Result is: 0.9432    0.7770    1.7070   11.8093
%e   = skt_cumulants(

if nargin<5;
    n=4;
end;

if nargin<4;
    df=Inf;
end;

if nargin<3;
    shape=0;
end;

if nargin<2;
    scale=1;
end;

if nargin<1;
    location=0;
end;

if isnan(n);
    n=4;
end;

if isnan(shape);
    shape=0;
end;

if isnan(scale);
    scale=1;
end;

if isnan(location);
    location=0;
end;

if isempty(n)
    n=4;
end;
if isempty(df)
    df=Inf;
end;
if isempty(shape)
    shape=0;
end;
if isempty(scale)
    scale=1;
end;
if isempty(location)
    location=0;
end;

%Define constraints on the scale parameter
if scale<=0
    error('The scale parameter must have a positive nonzero value');
end;
%Define constraints on the scalar value of n
if n<1 ||n>4||n>=df||round(n)~=n
    error('The value of n must be from 1 to 4 and smaller than df');
end;


if df==Inf;
    cum=sn_cumulants(shape,n);
elseif df==n;
    display('error: we need that df>n');
elseif df>=n;
    par=[location scale shape];
    delta=par(:,3)./sqrt(1+par(:,3).^2);
    mu=delta.*sqrt(df./pi).*exp(gammaln((df-1)./2)-gammaln(df./2));
    a=size(par);
    b=a(1,1);
    cum=zeros(b,n);
    cum(:,:)=NaN;
    cum(:,1)=mu;
    if n>1;
        cum(:,2)=df./(df-2)-mu.^2;
    end;
    if n>2;
        cum(:,3)=mu.*(df.*(3-delta.^2)/(df-3)-3.*df./(df-2)+2.*mu.^2);
    end;
    if n>3;
        cum(:,4)=(3.*df.^2/((df-2)*(df-4))-4.*mu.^2*df*(3-delta.^2)/(df-3)+6.*mu.^2.*df./(df-2)-3*mu.^4)-3*cum(:,2).^2;
    end;
    [app2,app1]=meshgrid(1:n,scale);
    kap=app1.^app2;
    %for i=1:n
    %   app(:,i)=(scale').^i;
    %end;
    cum=cum.*kap;
    cum(:,1)=cum(:,1)+par(:,1);
end;