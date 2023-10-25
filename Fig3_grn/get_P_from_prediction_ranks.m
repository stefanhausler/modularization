function [P] = get_P_from_prediction_ranks(prediction,pidx);
% GET_P_FROM_PREDICTION_RANKS  Significance levels of interactions according to ranks. 
%    P = GET_P_FROM_PREDICTION_RANKS(PREDICTION,PIDX) for a vector of ranks PIDX and interactions PREDICTION (#interactions-by-3) 
%    with rows [TF_INDEX, NTF_INDEX, PVALUE] 
%
%       TF_INDEX           index of regulating transcription factor (TF)
%       NTF_INDEX          index of regulated target gene (NTF)
%       PVALUE             p-value for rejecting the null hypothesis 'no interaction'
%
%    returns a vector P containing the p-values PVALUE of the interactions PREDICTION with ranks specified in PIDX.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

pidx(pidx>size(prediction,1)) = [];
pvalue = [0:1/size(prediction,1):1-1/size(prediction,1)];
P = pvalue(pidx); 

