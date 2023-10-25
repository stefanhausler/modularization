function [MAMA,MHMA,MH2MA] = precompute_moments(x,preNodes,postNodes,par);
% PRECOMPUTE_MOMENTS  Computes moment tensors used in CORRECT_GRN_NETWORK.
%    [MAMA,MHMA,MH2MA] = PRECOMPUTE_MOMENTS(X,PRENODES,POSTNODES,PAR) for gene expression levels X, ... and parameters PAR returns
%
%       MAMA,MHMA,MH2MA   moment tensors as used in CORRECT_GRN_NETWORK of size
%                         #(regulated target genes)-by-#(regulated target genes)-by-#(regulating transcription factors).  
%
%    X contains gene expression levels (#genes-by-#samples).
%
%    PRENODES is a vector (#(regulating transcription factors)-by-1) with the indices of all TF.
%
%    POSTNODES is a vector (#(regulated target genes)-by-1) with the indices of all NTF.
%
%    PAR contains parameters for the correction.
%
%       PAR.thresh        threshold parameter (1-by-1) for the test
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)




thresh = par.thresh;
nSamples = size(x,2);
nPreNodes = length(preNodes);
nPostNodes = length(postNodes);


tic
for i = 1:nPreNodes
    count = i;
    COUNT = nPreNodes;
    fprintf('Pre compute moments for moment ratios: %g%%, %g s         \r',100*count/COUNT,toc/(count-1)*(COUNT-count+1))
    
    y = x(preNodes(i)*ones(1,nPostNodes),:).*x(postNodes,:); % correlation of each prenode signale and all postnode signals
    
    lccmean = mean(y,2);
    lccstd = std(y,[],2);
    rho = corrcoef(y');
    
    % Notation as in Marsaglia 2006
    % z/w = 1/r * (a+x)/(b+y) + s
    
    mu_z = lccmean*ones(1,length(lccmean)); % all combinations of postnode pairs: mean and variance [me,v]
    mu_w = mu_z';
    
    sig_z = lccstd*ones(1,length(lccstd)) / sqrt(nSamples); % SEM
    sig_w = sig_z';
    
    % for correct limit calculation below
    rhoidx = find(abs(rho)>0.99999999);
    rho(rhoidx) = sign(rho(rhoidx))*0.99999999;
    
    b = mu_w./sig_w;
    a = (mu_z./sig_z-rho.*mu_w./sig_w)./sqrt(1-rho.^2);
    
    sign_r_a = (-1).^double((b >= 0) ~= (a >= 0));
    a = a.*sign_r_a;
    r = sig_w./(sig_z.*sqrt(1-rho.^2));
    r = r.*sign_r_a;
    s = rho.*sig_z./sig_w;
    
    % medians
    me = 1./r .* (a./b) + s;  % for correlated z and w. BUT: The median m and m2 are identical, only the variance changes.
    
    % variance
    % (2nd order approx. of the variance of the ratio (a+x)/(b+y) for x,y standard normals (see Kempen and van Vliet, 2000))
    v = (1./(b.^2) + a.^2./(b.^4));
    v = v ./ (r.^2); %  For correlated z and w: transform back to real distribution
    
    idx = find(abs(b) < thresh);
    v(idx) = Inf; % 0 contribution for inference below
    
    mtmp = 1./v;
    mtmp(1:1+size(mtmp,1):end) = 0; % remove diagonal elements
    MAMA(:,:,i) = mtmp;
    
    mtmp = me./v;
    mtmp(1:1+size(mtmp,1):end) = 0; % remove diagonal elements
    MHMA(:,:,i) = mtmp;
    
    mtmp = MHMA(:,:,i).*me;
    mtmp(1:1+size(mtmp,1):end) = 0; % remove diagonal elements
    MH2MA(:,:,i) = mtmp;
    
end % i
fprintf('\n');



