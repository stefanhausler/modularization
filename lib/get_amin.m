function [amin] = get_amin(par,nSites)
% GET_AMIN  Minimum significance level.
%    AMIN = GET_AMIN(PAR,NSITES) returns the minimum significance level in Equ. 19
%    as a function of the treshold PAR.THRESH and number of observable components
%    NSITES without the reference variable s_ref.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

nb = nSites*(nSites-1)/2;                                     % total number of ration moments
amin = 1-erf((par.thresh-1)/sqrt(2))^nb*erf(6/sqrt(2))^nb;    % Equ. 19
