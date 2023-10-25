function [pcc] = get_pcc_from_prediction_ranks(prediction,grnnet);
% GET_PCC_FROM_PREDICTION_RANKS Set p-values of interactions according to the ranks of evidence.
%    PCC = GET_PCC_FROM_PREDICTION_RANKS(PREDICTION,GRNNET) for gene regulatory network GRNNET and interactions PREDICTION
%    (#interactions-by-3) with rows [TF_INDEX, NTF_INDEX, EVIDENCE] sorted by EVIDENCE in descending order.
%
%       TF_INDEX           index of regulating transcription factor (TF)
%       NTF_INDEX          index of regulated target gene (NTF)
%       EVIDENCE           evidence for interaction
%
%       GRNNET.x           gene expression levels (#genes-by-#samples)
%
%    returns pseudo p-values PCC (#genes-by-#genes) for rejecting the null hypothesis 'no interaction'.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

[~,idx] = sort(-prediction(:,3));
if ~isequal(idx,[1:length(idx)]')
    warning('PREDICTION isn''t sorted!')
end

pvalue = [0:1/size(prediction,1):1-1/size(prediction,1)];
pcc = ones(size(grnnet.x,1));
pcc(prediction(:,2) + (prediction(:,1)-1)*size(pcc,2)) = pvalue;


