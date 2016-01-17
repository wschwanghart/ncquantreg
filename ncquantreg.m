function [p,h] = ncquantreg(x,y,varargin)

% Non-crossing polynomial quantile regression
%
% Syntax
%
%     p = ncquantreg(x,y)
%     p = ncquantreg(x,y,n,tau)
%     p = ncquantreg(x,y,n,tau,pn,pv,...)
%     [p,h] = ... 'plot',1 or 2)
%
% Description
%
%     ncquantreg finds the coefficients of a polynomial p(x) of degree n 
%     that fits the data in vector x to the quantiles tau of y. 
%     
%     ncquantreg(x,y) performs median regression (tau = 0.5) using a 
%     polynomial of degree n=1. 
%
%     ncquantreg(x,y,n,tau) fits a numel(tau) polynomials with degree n. 
%     The algorithm uses a stepwise multiple quantile regression estimation
%     using non-crossing constraints (Wu and Liu, 2009). The approach is
%     stepwise in a sense that a quantile function is estimated so that it
%     does not cross with a function fitted in a previous step. The
%     algorithm starts from the middle quantile (i.e. the one closest to
%     0.5) and than progressivly works through the quantiles with 
%     increasing distance from the middle.
%
%     ncquantreg(x,y,n,tau,pn,pv,...) takes several parameter name value
%     pairs that control the algorithm and plotting. 
% 
% Input arguments
%
%     x      independent variable (vector)
%     y      dependent variable
%     n      degree of polynomial (n = 1 (straight line, default)
%     tau    vector of quantiles (e.g. [0.1:0.1:0.2]). Values must range
%            between >0 and <1 (default tau = 0.5).
%
% Parameter name/value pairs
%
%     'xrange'   default [min(x) max(x)]. 
%                mx2 vector indicating the range of x that must fulfill
%                the constraint of non-crossing
%     'plot'     default 0
%                0: no plot
%                1: points in x and y and lines
%                2: lines only
%                
% Output
% 
%     p      n-by-numel(tau) matrix with polynomial coefficients.
%     h      handles to lineseries in plot if 'plot' is set to 1 or 2.
%            If 'plot' is 1, h(1) is the handle to the points.
%            h(2:numel(tau)+1) are the handles to the quantile curves.
%            If 'plot' is 2, h(1:numel(tau)) are the handles to the curves.
%
% Example
%
%     x = rand(100,1);
%     y = rand(100,1)+x*0.2+x.^2*.2;
%     n = 1;
%     tau = [0.1:0.1:0.9];
%     [b,h] = ncquantreg(x,y,n,tau,'xrange',[-1 2],'plot',1);
%     set(h(6),'LineWidth',2)
%
% Reference
% 
%     Wu, Y., Liu, Y., 2009. Stepwise multiple quantile regression
%     estimation using non-crossing constraints. Statistics and its
%     Interface 2, 299–310.
%
%     
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 17. January, 2016


% find minimum and maximum values
xmin = min(x);
xmax = max(x);

%% Parse Inputs
pp = inputParser;         
pp.FunctionName = 'ncquantreg';
addRequired(pp,'x');
addRequired(pp,'y',@(y) isequal(size(x),size(y)));
addOptional(pp,'n',1,@(n) validateattributes(n,{'numeric'},...
    {'scalar','>=',0,'integer','real','finite'},...
    'ncquantreg','n',3));
addOptional(pp,'tau',0.5,@(tau) validateattributes(tau,{'numeric'},...
    {'vector','>',0,'<',1,'real','finite'},...
    'ncquantreg','tau',4));

addParamValue(pp,'xrange',[xmin xmax],@(x) numel(x)==2 && x(1)<x(2));
addParamValue(pp,'plot',0,@(x) isscalar(x));

parse(pp,x,y,varargin{:});

%% Prepare
x = double(x(:));
y = double(y(:));
tau = pp.Results.tau;
n   = pp.Results.n;
ntau = numel(tau);

% design matrix
X = bsxfun(@power,x,0:n);

if isequal([xmin xmax],pp.Results.xrange);
    XCON = X;
else
    XCON = [X; bsxfun(@power,pp.Results.xrange',0:n)];
    XCON(x<pp.Results.xrange(1) | x>pp.Results.xrange(2),:) = [];    
%     XCON = linspace(xmin,xmax,100)';
%     XCON = bsxfun(@power,XCON,0:n);
end

% initial values guess is a least squares fit
p0 = X\y;

%% Quantile regression
if ntau == 1
    % if number of tau is one
    rho   = @(r)sum(abs(r.*(tau-(r<0)))); 
    p     = fminsearch(@(p)rho(y-X*p),p0);
else
    % if there is more than one tau
    tau = tau(:);
    % sort tau in ascending order
    [tau,ixf] = sort(tau,'ascend');   
    % find order to process the taus
    [~,ix]    = sort(abs(tau-0.5),'ascend');
    si        = sign((1:ntau)'-ix(1));
    % index into previous tau
    ixp       = ((1:ntau)'-si);
    
    % preallocate polynomials
    p     = zeros(n+1,ntau);
    options=optimset('Algorithm','active-set','Display','off');
    
    % go through all taus
    for r = 1:ntau;
        
        % anonymous function for residuals
        rho   = @(res)sum(abs(res.*(tau(ix(r))-(res<0))));
        
        if r == 1;
            % the first evaluation is unconstrained
            p(:,ix(r))     = fminsearch(@(p)rho(y-X*p),p0);
        else
            % all subsequent tau are constrained
            s = sign(tau(ix(r))-tau(ixp(ix(r))));
            % update initial values to coefficients of previous quantile
            p0         = p(:,ixp(ix(r)));
            p(:,ix(r)) = fmincon(@(p)rho(y-X*p),p0,...
                -s*XCON,-s*XCON*p0,...
                [],[],[],[],[],options);
        end
    end
    
    % set output coefficient matrix to user-supplied order of tau 
    p(:,ixf) = p;
end

%% Plot
if pp.Results.plot>0
    if ishold
        holdison = true;
    else
        holdison = false;
    end
    
    if pp.Results.plot == 1
        h = plot(x,y,'.');
        hold on
    else
        h = [];
    end
    
    xhat = linspace(min(x),max(x),100)';
    Xhat = bsxfun(@power,xhat,0:n);
    
    Yhat = Xhat*p;
    h = [h; plot(xhat,Yhat,'Color',[.6 .6 .6])];

    if ~holdison
        hold off
    end
else
    h = [];
end

