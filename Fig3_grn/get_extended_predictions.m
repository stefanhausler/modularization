function [ext_prediction] = get_extended_predictions(prediction,grnnet,varargin);
% GET_EXTENDED_PREDICTIONS  Extend prediction of gene interactions to 100000 interactions if less.
%    EXTENDED_PREDICTION = GET_EXTENDED_PREDICTIONS(PREDICTION,GRNNET) for gene regulatory network GRNNET and interactions PREDICTION
%    (#interactions-by-3) with rows [TF_INDEX, NTF_INDEX, EVIDENCE] sorted by EVIDENCE in descending order.
%
%       TF_INDEX           index of regulating transcription factor (TF)
%       NTF_INDEX          index of regulated target gene (NTF)
%       EVIDENCE           evidence for the interaction
%
%       GRNNET.x           gene expression levels (#genes-by-#samples)
%       GRNNET.TFUsed      indices of all TF (#(regulating transcription factors)-by-1)
%       GRNNET.NTFUsed     indices of all NTF (#(regulated target genes)-by-1)
%
%    returns EXTENDED_PREDICTION with the same format as PREDICTION but added interactions if less than 100000.
%
%    EXTENDED_PREDICTION = GET_EXTENDED_PREDICTIONS(PREDICTION,GRNNET,'only TF-NTF') removes in addition all interactions
%    that are not TF-NTF.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

preNodes = grnnet.TFUsed;
postNodes = grnnet.NTFUsed;
pre = preNodes( repmat([1:length(preNodes)],[length(postNodes) 1]) );
post = postNodes( repmat([1:length(postNodes)]',[1, length(preNodes)]) );

if (nargin > 2) & isequal(varargin{1},'only TF-NTF')
    
    % to keep order of submitted list which might contain identical values
    pcc = ones(size(grnnet.x,1));
    pcc(prediction(:,2) + (prediction(:,1)-1)*size(pcc,2)) = [0:1/size(prediction,1):1-1/size(prediction,1)];
    p = pcc(postNodes,preNodes);
    
    % to keep original value
    pcc2 = ones(size(grnnet.x,1));
    pcc2(prediction(:,2) + (prediction(:,1)-1)*size(pcc,2)) =  1-prediction(:,3);
    p2 = pcc2(postNodes,preNodes);
    
    
    [~, idx] = sort(p(:));
    ext_prediction = [pre(idx), post(idx), 1-p2(idx)];
    
else
    
    % all possible interactions sorted
    pcc = zeros(size(grnnet.x,1));
    p = pcc(postNodes,preNodes);
    [ps, idx] = sort(p(:));
    gsidx = [pre(idx), post(idx)];
    [ext_prediction,idx2] = unique([prediction(:,1:2); gsidx],'rows','stable');
    ps = [prediction(:,3); zeros(size(gsidx,1),1)];
    ext_prediction(:,3) = ps(idx2);
    
end

ext_prediction = truncate_predictions(ext_prediction);


