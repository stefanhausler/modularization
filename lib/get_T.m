function [p,T,zeta] = get_T(stat,idx,par)
% GET_T  Test statistic for modularization.
%    [P,T,ZETA] = GET_T(STAT,IDX,PAR) for moment ratio matrix and covariance matrix in STAT (obtained from GET_B), concise
%    index sets IDX (obtained from GET_IDX) and correction par.METHOD return
%
%       P             p-pvalue for the modularization corresponding to the concise index sets IDX
%       T             test statistic for the modularization corresponding to the concise index sets IDX
%       ZETA          shape parameter of the associated gamma distribution
%
%    PAR.METHOD determines the implemented correction.
%
%       METHOD == 0   no correction
%       METHOD == 1   simple correction
%       METHOD == 2   advanced correction
%
%    The significance levels for P have to be corrected according to Equ. 16 and 17 in the article.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

       METHOD = par.METHOD;

       if ~any(METHOD==[0 1 2])
          error('METHOD must be 0, 1 or 2!')
       end
       
       nH = size(idx,2); 
       
       bH = stat.b(max(idx,1)').*(idx'~=0);;
       sig2H = stat.sig2(max(idx,1)').*(idx'~=0);
       idx0H = (stat.idx0(max(idx,1)').*(idx'~=0))==1;

       if METHOD==1
         i = find(idx(:)); % take only non-zero indices
         SIG0 = stat.SIG(i,i); 
         SIGd = diag(diag(1./sqrt(SIG0)));
         [~,LAMMAX] = eig( SIGd*SIG0*SIGd);
         LAMMAX = max([LAMMAX(:); 1]); % to prevent numeric problems
         sig2H = LAMMAX*sig2H;

       elseif METHOD==2

         m = [];
         SIG0T = [];
         for i5 = 1:nH
            i = find(idx(:,i5)); % skip zero indices, which allows different row lengths
            SIG0p = stat.SIG(idx(i,i5),idx(i,i5)); % take only non-zero indices

            [Q,lam] = eig(SIG0p);
            rr = sum(Q,1)';
            sig2H(i5,i) = (rr.^-2).*diag(lam); % only for non-zero indices
            bH(i5,i) = diag(rr.^-1)*Q'*bH(i5,i)'; % only for non-zero indices

            m = blkdiag(m,lam^(-1/2)*Q');
         end
         i = find(idx(:));
         SIG0T = stat.SIG(idx(i),idx(i)); % take only non-zero indices
         
         M = m*SIG0T*m';
         [~,lamM] = eig((M+M')/2);   % to prevent numeric imag problems due to machine precision
         LAMMAX = max([lamM(:); 1]); % 1 needed if M is []
         
         sig2H = LAMMAX*sig2H;
       end

 
       % clamp b to 0
       bH(idx0H) = 0;
       sig2H(idx0H) = Inf;
       sig2H(idx'==0) = Inf; % don't sum over zero indices below

       % no clamping of singletons because they dont't enter the statistic anyway
       den = sum(1./sig2H,2);
       te = sum( bH.^2./sig2H,2) - sum( (bH./sig2H),2).^2 ./ den;  
   
       % ensure correct asymptotics
       te(den==0) = 0; 

       zeta = 1/2*sum(sum(idx~=0) - 1);  % shape parameter of the associated gamma distribution
       T = 1/2*sum(te);                  % test statistic
       logit = T - zeta;
       
       p = gammainc(max(T,0),zeta,'upper'); % max to prevent numeric problems

