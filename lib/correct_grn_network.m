function [gsidx,p] = correct_grn_network(grnnet,par);
% CORRECT_GRN_NETWORK  Correction to the evidence for gene interactions predicted from gene expression data.
%    [GSIDX,P] = CORRECT_GRN_NETWORK(GRNNET,PAR) for gene expression data GRNNET and correction parameters PAR returns
%
%       GSIDX             predictions of interactions (#interactions-by-3) with rows 
%                         [TF_INDEX, NTF_INDEX, EVIDENCE] sorted by EVIDENCE in descending order.
%
%                         TF_INDEX           index of regulating transcription factor (TF)
%                         NTF_INDEX          index of regulated target gene (NTF)
%                         EVIDENCE           evidence for the interaction (0>=EVIDENCE>=1)
%                                            If GRNNET.PCC are p-values, then (1-EVIDENCE) are p-values
%                                            for rejecting the null hypothesis 'no interaction'.
%
%       P                 EVIDENCE (#(regulated target genes)-by-#(regulating transcription factors)) for all 
%                         TF-NTF interactions
%
%    GRNNET contains gene expression data.
%
%       GRNNET.x          gene expression levels (#genes-by-#samples)
%       GRNNET.TFUsed     indices of all TF (#(regulating transcription factors)-by-1)
%       GRNNET.NTFUsed    indices of all NTF (#(regulated target genes)-by-1)
%
%    PAR contains uncorrected evidence of interactions and parameters for the correction.
%
%       PAR.thresh        threshold parameter (1-by-1) for the test (default 0)
%       PAR.pcc           evidence (#genes-by-#genes) for no interactions (0>=PAR.pcc>=1),
%                         e.g. p-values for rejecting the null hypothesis 'no interaction' 
%                         or rescaled interaction ranks (equal values allowed)
%       PAR.alpha         significance level (1-by-1) for PAR.pcc to select interactions to be tested for correction
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)



% verify existence of MEX file minidx.c
if exist('minidx')~=3
    disp('Compiling MEX file minidx.c.')
    mex('-setup','C')
    st = dbstack;
    libpath = fileparts(which(st(1).name));
    cdir = pwd;
    try
        chdir(libpath)
        mex('minidx.c')
        chdir(cdir)
    catch
        chdir(cdir)
        error('Problem compiling MEX file minidx.c!')
    end
end

% initialization
thresh = par.thresh;
x = grnnet.x;
nGenes= size(x,1);
nSamples = size(x,2);
pcc = par.pcc;
alpha = par.alpha;
preNodes = grnnet.TFUsed;
postNodes = grnnet.NTFUsed;
nPreNodes = length(preNodes);
nPostNodes = length(postNodes);

% for conversion from nGenes x nGenes format to nPostNodes x nPreNodes index
cPreNodes = zeros(max(preNodes),1);
cPostNodes = zeros(max(postNodes),1);
cPreNodes(preNodes) = 1:length(preNodes);
cPostNodes(postNodes) = 1:length(postNodes);


% find edges included in the test
pW = (pcc < alpha);

nSamplesh = floor(nSamples/2);

fit_setup_tmp.x = sum(x);
fit_setup_tmp.size_x = size(x);
fit_setup_tmp.nSamplesh = nSamplesh;
fit_setup_tmp.preNodes = preNodes;
fit_setup_tmp.postNodes = postNodes;
fit_setup_tmp.thresh = par.thresh;

persistent fit_setup MAMA MHMA MH2MA MAMAb MHMAb MH2MAb

if isempty(fit_setup)|~isequal(fit_setup,fit_setup_tmp)
    fit_setup = fit_setup_tmp;
    
    MAMA = []; MHMA = []; MH2MA = []; MAMAb = []; MHMAb = []; MH2MAb = [];  % in case variables are already set
    
    [MAMA,MHMA,MH2MA] = precompute_moments(x(:,1:2:2*nSamplesh),preNodes,postNodes,par);
    [MAMAb,MHMAb,MH2MAb] = precompute_moments(x(:,2:2:2*nSamplesh),preNodes,postNodes,par);
end


nans1 = sum(isnan(MAMA(:)));
nans2 = sum(isnan(MHMA(:)));
nans3 = sum(isnan(MH2MA(:)));

if (nans1+nans2+nans3)>0
    fprintf('Warning: %g feature matrix elements are not a number!\n',max([nans1,nans2,nans3])/759)
end


resN = zeros(nPostNodes,nPreNodes);
resL = Inf*ones(nPostNodes,nPreNodes);


% iterate over transcription factors
for i = 1:nPreNodes
    
    % find all valid post nodes
    m = postNodes( find(pW(postNodes,preNodes(i))) );
    m(m<=preNodes(i)) = []; % c) only links to nodes with larger indices if pre and post nodes are overlapping
    
    % sort with most significant links first (link 1), only check motifs with less significant links 2,3,4.
    [~,k] = sort( pcc(m,preNodes(i)) );
    m = m(k);
    
    if size(m,2) ~= 1
        m = m';
    end
    
    % potential post nodes connected to preNodes(i), i.e. link 1 and 4
    pn1 = m(:,ones(1,length(m)));
    pn2 = pn1';
    pn1 = pn1(:);
    pn2 = pn2(:);
    j = find(pn2>pn1); % because pn1 has lower p-value it has lower index
    pn1 = pn1(j);
    pn2 = pn2(j);
    
    if isempty(pn1)
        continue
    end
    
    hc4 = pcc(pn1,preNodes(i)); % p-value of link 1
    
    % loop through all preNodes with
    %  a) except i: check also prenodes smaller than i, because the preNodes are not sorted according to correltion strength
    %  b) connected with at least one link to preNodes(j) (and both to preNodes(i))
    %  c) both post-nodes indices larger then both pre-nodes indices (to keep directionality of links)
    %  d) correlations of links 2 and 3 are smaller than for link 1
    
    for j = setdiff(1:nPreNodes,i) % a)
        
        fprintf('Inference: %g/%g %g/%g              \r',i,nPreNodes,j,nPreNodes-1)
        
        % search for links between prenode(j) and those postnodes connected to prenode(i) -> three links exist
        c1 = pW(pn1,preNodes(j)) | pW(pn2,preNodes(j)); % b)
        c2 = (preNodes(j) < pn1) & (preNodes(j) < pn2); % c)
        c3 = (pcc(pn1,preNodes(j)) > hc4) & (pcc(pn2,preNodes(j)) > hc4);   % d) ... > because p-values are opposite of correlations
        % c = find(c1 & c2 & c3);
        c4 = (pcc(pn1,preNodes(j)) >= hc4) & (pcc(pn2,preNodes(j)) >= hc4) & (i<j); % d) ... equal p-valued links to preNodes(j) allowed if i < j to
        % still check e.g. all p=0 motifs but drop ties
        c = find( c1 & c2 & (c3|c4) );
        
        % if no testable motifs are found then skip to next iteration
        if isempty(c)
            continue
        end
        
        % convert true node index back to the corresponding index in the smaller nPostNodes x nPreNodes format
        pn1c = cPostNodes(pn1(c));
        pn2c = cPostNodes(pn2(c));
        
        % check motif for nodes preNodes(i) preNodes(j) pn1 pn2;
        
        maxLogitDiff = zeros(length(pn1c),1);
        
        % check 1/4 = 3/2
        hpn1 = pn1c + (pn2c-1)*nPostNodes;
        i1 = hpn1 + (i-1)*nPostNodes*nPostNodes; % indices for link 1/link 4 ratio
        i2 = hpn1 + (j-1)*nPostNodes*nPostNodes; % indices for link 3/link 2 ratio
        
        % check 4/1 = 2/3
        hpn2 = pn2c + (pn1c-1)*nPostNodes;
        i3 = hpn2 + (i-1)*nPostNodes*nPostNodes; % indices for link 4/link 1 ratio
        i4 = hpn2 + (j-1)*nPostNodes*nPostNodes; % indices for link 2/link 3 ratio
        
        mAmA1  = [MAMA(i1) MAMAb(i2)];
        mHmA1  = [MHMA(i1) MHMAb(i2)];
        mH2mA1 = [MH2MA(i1) MH2MAb(i2)];
        
        mAmA2  = [MAMA(i3) MAMAb(i4)];
        mHmA2  = [MHMA(i3) MHMAb(i4)];
        mH2mA2 = [MH2MA(i3) MH2MAb(i4)];
        
        
        %%%%%%%%%%%%%%%%%
        % compute logits
        
        logitAllNoise = sum(mHmA1,2).^2 ./ sum(mAmA1,2) / 2;
        logitAllStructure = sum ( mH2mA1,2 ) / 2;
        logitDiff = logitAllStructure - logitAllNoise;
        logitDiff(isnan(logitDiff)) = Inf;  % take out of analysis. Results is p-pvalue of 0 -> no impact on the selection based on cc
        logitDiff1 = max(logitDiff,0);      % must be > 0 arithmetically (see proof), below => errors due to machine precision
        
        
        logitAllNoise = sum(mHmA2,2).^2 ./ sum(mAmA2,2) / 2;
        logitAllStructure = sum ( mH2mA2,2 ) / 2;
        logitDiff = logitAllStructure - logitAllNoise;
        logitDiff(isnan(logitDiff)) = Inf;  % take out of analysis. Results is p-pvalue of 0 -> no impact on the selection based on cc
        logitDiff2 = max(logitDiff,0);      % must be > 0 arithmetically (see proof), below => errors due to machine precision
        
        [maxLogitDiff] = max([logitDiff1 logitDiff2 maxLogitDiff],[],2);
        
        % search for least significant result for motif tests including links 2 and 3 (from pn1 and pn2 to preNodes(j)), respectively
        [rl,rn] = minidx(maxLogitDiff,pn1c,pn2c,nPostNodes);
        rl(rl==1e20) = Inf;
        
        % store least significant result for motif tests including links 2 and 3 (from pn1 and pn2 to preNodes(j))
        resL(:,j) = min([resL(:,j) rl],[],2);
        resN(:,j) = resN(:,j) + rn;
        
    end % j
    
end % i

pF = 1 - gammainc(resL,0.5);

% destroy potential order of equivalent pF == 1
i = find(pF==1);
pF(i) = pF(i) - 1e-5*rand(size(pF(i)));

pL = pcc(postNodes,preNodes);

p = max(pL,alpha*pF);

% extract list of TF-NTF interactions sorted by evidence from p

[ps, idx] = sort(p(:));

pre = preNodes( repmat([1:nPreNodes],[nPostNodes 1]) );
post = postNodes( repmat([1:nPostNodes]',[1, nPreNodes]) );

gsidx = [pre(idx), post(idx), 1-ps];
i = find( ~ismember(gsidx(:,1),preNodes) | ~ismember(gsidx(:,2),postNodes));
gsidx(i,:) = [];








