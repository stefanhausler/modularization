function [gsidx] = get_prediction_from_pcc(pcc,grnnet);
% GET_PREDICTION_FROM_PCC  Sort TF-NTF gene interactions according to their Pearson correlation coefficients.
%    GSIDX = GET_PREDICTION_FROM_PCC(PCC,GRNNET) returns interactions GSIDX (#interactions-by-3) with rows
%    [TF_INDEX, NTF_INDEX, PVALUE] sorted by the PVALUE in ascending order. 
%
%       TF_INDEX           index of regulating transcription factor (TF)
%       NTF_INDEX          index of regulated target gene (NTF)
%       PVALUE             p-value for rejecting the null hypothesis 'no interaction'
%
%    PCC p-values (#genes-by-#genes) for rejecting the null hypothesis 'no interaction'.
%
%    GRNNET contains gene expression data.
%
%       GRNNET.TFUsed     indices of all TF (#(regulating transcription factors)-by-1)
%       GRNNET.NTFUsed    indices of all NTF (#(regulated target genes)-by-1)
%
%    Contraints for extracted TF-NTF gene interactions are i) no self interactions and ii)  only TF can regulate.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

preNodes = grnnet.TFUsed;
postNodes = grnnet.NTFUsed;
nPreNodes = length(preNodes);
nPostNodes = length(postNodes);

pcc = pcc + eye(size(pcc)); % no self interactions
pcc(:,nPreNodes+1:end) = 1; % only TF can regulate
[ps, idx] = sort(pcc(:));

pre = repmat([1:size(pcc,2)],[size(pcc,1) 1]);
post = repmat([1:size(pcc,1)]',[1, size(pcc,2)]);

gsidx = [pre(idx), post(idx), ps];


