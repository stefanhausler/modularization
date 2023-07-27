function [idx] = get_idx_linear(S,nSites)
% GET_IDX_LINEAR  Concise index set for a single linear functional module.
%    IDX = GET_IDX_LINEAR(S,NSITES) for functional module S and total number of observable components NSITES
%    (without the reference variable s_ref) returns the indices of the elements of the NSITES-by-NSITES moment
%    ratio matrix (obtained by GET_B) that correspond to the concise index sets for the compoments within the
%    functional modults S if these contribute only linearly to the overall functional organization.
%
%    IDX is an input argument for the function GET_T.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

       i = zeros(nSites);
       i(S,S) = 1;
       i = triu(i ~= 0,1);
       idx = find(i);
