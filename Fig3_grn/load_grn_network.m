function [grnnet] = load_grn_network(path_dream5,NETWORK,varargin);
% LOAD_GRN_NETWORK  Load gene regulatory network of the DREAM5 challenge.
%    GRNNET = LOAD_GRN_NETWORK(PATH_DREAM5,NETWORK) for path PATH_DREAM5 and index NETWORK returns
%
%       GRNNET.x          gene expression levels (#genes-by-#samples)
%       GRNNET.TFUsed     indices of all TF (#(regulating transcription factors)-by-1)
%       GRNNET.NTFUsed    indices of all NTF (#(regulated target genes)-by-1)
%       GRNNET.gsidx      gold standard (#interactions-by-3)
%
%    PATH_DREAM5 (string) specifies the location of the DREAM5 folder.
%
%    NETWORK index 1 ... In silico network
%                  3 ... E. coli network
%
%    GRNNET = LOAD_GRN_NETWORK(PATH_DREAM5,NETWORK,FOLDS) loads a subsample of the gene expression levels
%    specified in the binary vector FOLDS. Non-zero Elements of FOLDS indicate selected samples. If the
%    length of FOLDS is smaller than the number of samples, FOLDS is expanded by repeated copying.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)


[x,info] = get_grn_data(NETWORK,path_dream5);

if (nargin == 2)
    fidx = 1:size(x,2);
else
    FOLDS = varargin{1};
    if size(FOLDS,1)>size(FOLDS,2)
        folds = FOLDS;
    elseif size(FOLDS,1)<size(FOLDS,2)
        folds = FOLDS';
    else
        error('Incorrect FOLDS.');
    end
    folds = repmat(folds,[ceil(size(x,2)/length(folds)) 1]);
    fidx = find(folds);
    fidx(fidx>size(x,2)) = [];
end

x = x(:,fidx);

% normalize data to zero mean and unit variance

nSamples = size(x,2);
x = x - repmat(mean(x,2),[1 nSamples]);
sd = std(x,[],2);
sd(isnan(sd)) = 1;
x = x./repmat(sd,[1 nSamples]);


nGenes = size(x,1);

% use only links with strong evidence as in the dream5 challenge
ridx = find(info.gsidx(:,3)==1);
gsidx = info.gsidx(ridx,1:2);

tfidx = info.tfidx(:,1);
ntfidx = setdiff([1:nGenes]',tfidx);

genesUsed = [1:nGenes]';
TFUsed = intersect(tfidx,genesUsed);
NTFUsed = setdiff(genesUsed,TFUsed);

gsidx = [gsidx, ones(size(gsidx,1),1)];
grnnet.x = x;
grnnet.TFUsed = TFUsed;
grnnet.NTFUsed = NTFUsed;
grnnet.gsidx = gsidx;



