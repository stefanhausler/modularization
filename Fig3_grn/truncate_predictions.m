function [gsidx] = truncate_predictions(gsidx);
% TRUNCATE_PREDICTIONS  Truncate gene interactions to at most 100000 interactions.
%    GSIDX = TRUNCATE_PREDICTIONS(GSIDX) for gene interactions GSIDX returns
%
%       GSIDX             predictions of interactions (#interactions-by-3) with rows 
%                         [TF_INDEX, NTF_INDEX, EVIDENCE]
%
%                         TF_INDEX           index of regulating transcription factor (TF)
%                         NTF_INDEX          index of regulated target gene (NTF)
%                         EVIDENCE           evidence for the interaction (0>=EVIDENCE>=1)
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)
           
gsidx(100001:end,:) = [];

