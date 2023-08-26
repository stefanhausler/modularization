% Statistical test for functional modularizations
%
% Example of how to use the toolbox for the statistical test in 
% "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% Runtime on an Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (28 cores):  < 1 sec execution time and < 1 GB memory
%
% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)


clear all

addpath('./lib')    % path to test library


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test parameters

par.thresh = 5;         % threshold of 5, 6 or 7
par.METHOD = 2;         % 0 ... uncorrected, 1 ... simple correction, 2 ... advanced correction

disp('Load data matrices.')

load('data.mat','x','z')     
% x ... (no. of network components without s_ref) -by- (no. of samples) 
%       samples of network states excl. the network component s_ref
% z ... 1 -by- (no. of samples)
%       samples of the network component s_ref

nSites = size(x,1);     % number of network components (without s_ref)
nS  = size(x,2);        % number of samples

disp('Normalize to mean zero and unit variance.')

x = x - repmat(mean(x,2),[1 nS]);
x = x./repmat(std(x,[],2),[1 nS]);
z = z - repmat(mean(z,2),[1 nS]);
z = z./repmat(std(z,[],2),[1 nS]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statistical test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Generate modularization and concise index sets.')

M = {[2,3,4],[3,4]};          % a nested modularization consisting of the two function modules [2,3,4] and [3,4]
IDX = get_idx(M,nSites);      % indices for concise index sets for the modularization M
IDXL = get_idx_linear(1:nSites,nSites);     % indices for concise index sets of M_L

% Evaluate modularizations

stat = get_B(z,x,par);    % get moment ratios, variances, covariance matrix and indices of elememts to clamp
p = get_T(stat,IDX,par);                % calculate the p-value for M
p_L = get_T(stat,IDXL,par);             % calculate the p-value for M_L

disp('Calculate Bonferroni corrected individual significance level.')

amin = get_amin(par,nSites); % minimum significance level
alpha = 0.01;                % nominal significance level
nH = 2;                      % number of hypotheses
a = (alpha-amin)/(1-amin);   % individual significance level
alpha_corr = a/nH;           % Bonferroni corrected individual significance level

fprintf('\nThe p-value for M is %g and the p-value for M_L is %g.',[p p_L]);
fprintf('\nThe Bonferroni corrected individual significance level is %g.\n',alpha_corr);
