function [stat] = get_B(sref,s,par)
% GET_B  Moment ratios to test for modularizations.
%    STAT = GET_B(SREF,S,PAR) for reference variable SREF (1-by-sample), observable components S (component-by-sample) and 
%    threshold PAR.THRESH returns
%
%       STAT.B       moment ratios matrix
%       STAT.SIG2    variances matrix
%       STAT.SIG     covariance matrix
%       STAT.IDX0    indices of elememts of B and SIG2 to be clamped because of the threshold
%
%    STAT is an input argument for the function GET_T.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

    thresh = par.thresh;  % threshold for clamping
    
    if ~any(thresh==[5 6 7])
       error('Threshold must be 5, 6 or 7!')
    end
    
    PARAM = [1.36686    2.0470    4.7350   -1.9230   -1.2310    2.7900;  
             1.38000    1.1760   18.9120  -19.2240    0.0220   26.6800;   
             1.53300    0.4460   39.5820  -65.2340   42.5320  120.6800];  

    PAR = PARAM(thresh-4,:); % correction parameters of function lambda (Equ. 15)
    

    nSites = size(s,1);      % number of sites

    % calculate numerators and denominators of the moment ratios
    t1 = 1:2:length(sref);
    u1 = s(:,t1);
    u2 = s(:,t1).*repmat(sref(1,t1),[nSites 1]);
    U = zeros(nSites,nSites,length(t1));
    for i1 = 1:length(t1)
       U(:,:,i1) = u1(:,i1)*u2(:,i1)';
    end
    
    t2 = 1:(length(sref)/2);
    t3 = (length(sref)/2+1):length(sref);
    w1 = s(:,t2).*repmat(sref(1,t2),[nSites 1]);
    w2 = s(:,t3).*repmat(sref(1,t3),[nSites 1]);
    W = zeros(nSites,nSites,length(t2));
    for i1 = 1:length(t2)
       W(:,:,i1) = w1(:,i1)*w2(:,i1)';
       W(:,:,i1) = triu(W(:,:,i1)) + triu(W(:,:,i1),1)'; % doesn't matter which of the two estimates (above or below the diagonal) is used
    end
  
    n = size(U,3);
   
    mu_u = mean(U,3);
    sig_u = std(U,[],3)/sqrt(n);
    mu_w = mean(W,3);
    sig_w = std(W,[],3)/sqrt(n);

    rho = mean((U-repmat(mu_u,[1 1 n])).*(W-repmat(mu_w,[1 1 n])),3)./(sig_u.*sig_w)/n; % n compensate 1/sqrt(n) in sig_u,sig_w to get real std

    % Calculate b and its variance in the standard form of ratios of normal distributions (Marsaglia et al., 2006)

    tb = mu_w./sig_w;
    ta = (mu_u./sig_u-rho.*mu_w./sig_w)./sqrt(1-rho.^2);
    sign_r_a = (-1).^double((tb >= 0) ~= (ta >= 0));
    ta = ta.*sign_r_a;
    r = sig_w./(sig_u.*sqrt(1-rho.^2));
    r = r.*sign_r_a;
    s = rho.*sig_u./sig_w;

    % variance
    vnoise2ndO = (1./(tb.^2) + ta.^2./(tb.^4)); % 2nd order approx. of the variance of the ratio (a+x)/(b+y) for x,y standard normals (see Kempen and van Vliet, 2000)

    % correction for the non-asymptotic test
    lambda = (1 + ( PAR(2)*abs(tb/6).^-2.0 + PAR(3)*abs(tb/6).^-4.0 + PAR(4)*abs(tb/6).^-6.0 + PAR(5)*abs(tb/6).^-8.0 + PAR(6)*abs(tb/6).^-10.0 ) )*PAR(1);

    sig2 = vnoise2ndO.* lambda;
    sig2 = sig2 ./ (r.^2); %  For correlated z and w: transform back to real distribution

    % actual distribution 
    b = mu_u./mu_w; 

    % transform to vector representations
    i1 = find(ones(nSites));
    
    Uv = permute(U,[3 1 2]);
    Wv = permute(W,[3 1 2]);
    Uv = Uv(:,i1);
    Wv = Wv(:,i1);
    Uv = permute(Uv,[2 1]);
    Wv = permute(Wv,[2 1]);
    nv = size(Uv,1);
    
    mu_uv = mu_u(i1);
    mu_wv = mu_w(i1);
    lambdav = lambda(i1);
    tbv = tb(i1);

    % full COV calculation (Equ. 11)
 
    SIGS = cov(cat(1,Uv,Wv)');
    SIGSS = SIGS;
    
    SIG11 = SIGS(1:nv,1:nv);
    SIG12 = SIGS(1:nv,nv+1:2*nv);
    SIG21 = SIGS(nv+1:2*nv,1:nv);
    SIG22 = SIGS(nv+1:2*nv,nv+1:2*nv);
   
    t1 = mu_wv*mu_wv';
    t2 = ones(nv,1)*(mu_uv'./mu_wv');
    t3 = (mu_uv./mu_wv)*ones(1,nv);
    t4 = (mu_uv*mu_uv') ./ t1;
   
    SIG = (1./t1).*(SIG11 - t2.*SIG12 - t3.*SIG21 + t4.*SIG22) / n ;
    SIG = sqrt(lambdav*lambdav').*SIG; 
 

    % clamping of sigma (before the eigendecomposition in get_T)
    idx0 = (abs(tb) < thresh); 
    idx0v = (abs(tbv) < thresh); 
    M = logical(repmat(idx0v,[1 size(idx0v,1)]) + repmat(idx0v',[size(idx0v,1) 1]));
    Md = logical(M.*eye(nv));

    SIG(M) = 0;   % set clamped off-diagonal elements to 0
    SIG(Md) = 1;  % set clamped on-diagonal elements to 1
    
    % output arguments
    stat.b = b;
    stat.sig2 = sig2;
    stat.idx0 = idx0;
    stat.SIG = SIG;
