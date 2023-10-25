% Statistical test for indirect interactions in gene regulatory networks
%
% Reproduces results shown in Fig 3 discussed in the section "Inferring functional modules in gene regulatory networks" in
% "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler
%
% Core function:
%
%   % gene regulatory network (all fields double)
%
%   grnnet.x:       [#genes x #samples]                        ... gene expression data
%   grnnet.TFUsed:  [#(regulating transcription factors) x 1]  ... indices of regulating transcription factors (TF)
%   grnnet.NTFUsed: [#(regulated target genes) x 1]            ... indices of regulated target genes (NTF)
%
%   % test parameters (all fields double)
%
%   par.thresh: [1 x 1]           ... threshold parameter of the test (default 0)
%   par.pcc:    [#genes x #genes] ... p-values for gene interactions and null hypothesis 'no interaction'
%   par.alpha:  [1 x 1]           ... significance level indicating interactions to consider for the correction
%
%   % apply method to correct par.pcc according to the statistical test
%
%   gsidx = correct_grn_network(grnnet,par);
%
%   % result
%
%   gsidx:  [#interactions x 3]      ... predictions with rows [TF index, NTF index, evidence], where 0>=evidence>=1

% Runtime on an Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (28 cores):  NETWORK=1;  < 15 min execution time and <  24 GB memory
%                                                                      NETWORK=3;  < 10 h   execution time and < 310 GB memory
%
% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)

%----------------------------------------------------------
% initialization
%----------------------------------------------------------

clear all
addpath('../lib')    % path to test library
rng(1);              % for repeatability


%----------------------------------------------------------
% DREAM5 network inference challenge

path_dream5 = '/home/haeusler/Projects/hierarchy/publication/code_v3/Fig3_grn/dream5';

NETWORK = 1 % In silico network
% NETWORK = 3 %  E. coli network

% DREAM5 matlab evaluation scripts
addpath(fullfile(path_dream5,'Evaluation scripts/matlab'))

% DREAM5 predictions for 36 inference methods
files = dir(sprintf('%s/Network_predictions/Network%g/*.txt',path_dream5,NETWORK));             % filenames of 36 methods
files = fullfile(path_dream5,'Network_predictions',sprintf('Network%g',NETWORK),{files.name});  % add full path to filenames

% load precomputed probability densities for various metrics
pdf_aupr  = load(sprintf('%s/Evaluation scripts/INPUT/probability_densities/Network%g_AUPR.mat',path_dream5,NETWORK));
pdf_auroc = load(sprintf('%s/Evaluation scripts/INPUT/probability_densities/Network%g_AUROC.mat',path_dream5,NETWORK));



%--------------------------------------------------------------------------------
% Example 1:
% Implementation of the statistical test given a list of interactions sorted by p-values.
% P-values obtained using Pearson's correlation coefficient.
% NOTE: Resulting list of interactions differs slightly from method 'Correlation 1' of the DREAM5 challenge because equal values are sorted differently
%--------------------------------------------------------------------------------

grnnet = load_grn_network(path_dream5,NETWORK);     % load grn network
[~, pcc] = corrcoef(grnnet.x');                     % p-values for all interactions and null hypothesis "no interaction"
prediction = get_prediction_from_pcc(pcc,grnnet);   % predicted interactions sorted by evidence

ext_prediction = truncate_predictions(prediction);  % truncate to 10^5 interactions as done for the DREAM5 challenge
                                                    % (required to use precomputed probability densities for various metrics)
gold_edges = grnnet.gsidx;
% Equivalent to original code (which might cause MATLAB:16587: GLib-CRITICAL warning):
% goldfile = sprintf('%s/Evaluation scripts/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network%g.tsv',path_dream5,NETWORK);
% gold_edges = load_dream_network(goldfile);


%------------------------------------------------------------------------------------------
% Evaluate the predictions (without correction by the statistical test)

[tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, ext_prediction, pdf_aupr, pdf_auroc);

perf_total = -(log10(p_auroc)+log10(p_aupr))/2;
fprintf('\nInference method Correlation 1 without correction: AUPR %g   AUROC %g   Total %g\n\n',aupr*100,auroc*100,perf_total)


%------------------------------------------------------------------------------------------
% Evaluate the predictions corrected by the statistical test

% set test parameters
par.thresh = 0;      % ... threshold parameter of the test
par.pcc = pcc;       % ... p-values for gene interactions and null hypothesis 'no interaction'
par.alpha = 10^-10;  % ... significance level indicating interactions to consider for the correction

% apply inference method
gsidx = correct_grn_network(grnnet,par);

% truncate predictions to 100000 interactions
gsidx = truncate_predictions(gsidx);

[tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, gsidx, pdf_aupr, pdf_auroc);

perf_total = -(log10(p_auroc)+log10(p_aupr))/2;
fprintf('\nInference method Correlation 1 with correction: AUPR %g   AUROC %g   Total %g\n\n',aupr*100,auroc*100,perf_total)



%----------------------------------------------------------------------------
% Example 2:
% Implementation of the statistical test given a list of ranked interactions.
% Interactions ranks obtained using predictions of DREAM5 inference method H0 (1...36)
% Results shown in Fig.3E-G
%----------------------------------------------------------------------------

% load grn network
grnnet = load_grn_network(path_dream5,NETWORK,[1 1 1 1 1 1 1 0]);         % don't select every 8th sample
grnnet_eval = load_grn_network(path_dream5,NETWORK,[0 0 0 0 0 0 0 1]);    % only select every 8th sample

gold_edges = grnnet.gsidx;                                                % load gold standard
% Equivalent to original code (which might cause MATLAB:16587: GLib-CRITICAL) warnings:
% goldfile = sprintf('%s/Evaluation scripts/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network%g.tsv',path_dream5,NETWORK);
% gold_edges = load_dream_network(goldfile);

midx = 21; % index of DREAM5_NetworkInference_Other1 (Genie3)
% midx = 8;  % index of DREAM5_NetworkInference_Correlation1
% midx = 7;  % index of DREAM5_NetworkInference_Community_Network

% load predictions of DREAM5 inference method
prediction = load_dream_network(files{midx});
pcc = get_pcc_from_prediction_ranks(prediction,grnnet);
P = get_P_from_prediction_ranks(prediction,[10:10:40 50:50:250 300:200:5000 6000:1000:10000 20000:10000:50000]);


%------------------------------------------------------------------------------------------
% evaluate predictions of the DREAM5 inference method midx using all predicted interactions

% extended predictions to 100000 interactions to avoid bias for shorter lists
% (induced by the probability densities for various metrics that have been precomputed for 100000 interactions)
ext_prediction = get_extended_predictions(prediction,grnnet);

[tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, ext_prediction, pdf_aupr, pdf_auroc);

perf_total = -(log10(p_auroc)+log10(p_aupr))/2;
fprintf('\nInference method %g using all predicted interactions: AUPR %g   AUROC %g   Total %g\n\n',midx,aupr*100,auroc*100,perf_total)


%--------------------------------------------------------------------------------------------------
% evaluate predictions of the DREAM5 inference method midx using only predicted TF-NTF interactions

% extended predictions to 100000 interactions to avoid bias for shorter lists and
% (induced by the probability densities for various metrics that have been precomputed for 100000 interactions)
ext_prediction = get_extended_predictions(prediction,grnnet,'only TF-NTF');

[tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, ext_prediction, pdf_aupr, pdf_auroc);

perf_total = -(log10(p_auroc)+log10(p_aupr))/2;
fprintf('\nInference method %g using only TF-NTF interactions without correction: AUPR %g   AUROC %g   Total %g\n\n',midx,aupr*100,auroc*100,perf_total)


%--------------------------------------------------------------------------------------------------
% evaluate predictions of the DREAM5 inference method midx on the hold out set for the gene expression levels

clear x AUPR
for nP = 1:length(P)
    fprintf('Network for correction %g/%g                           \n',nP,length(P))
    
    % set test parameters
    par.thresh = 0;
    par.alpha = P(nP);
    par.pcc = pcc;          % the same for grnnet_eval and grnnet
    
    % apply inference method
    gsidx = correct_grn_network(grnnet_eval,par);
    
    % truncate predictions to 100000 interactions
    gsidx = truncate_predictions(gsidx);
    
    % evaluate prediction using various metrics
    [tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, gsidx, pdf_aupr, pdf_auroc);
    
    AUPR(nP) = aupr;
end


%-----------------------------------------------------------------------------------------------------------------------
% evaluate predictions of the DREAM5 inference method midx on the gene expression levels not included in the hold out set

[~,nP] = max(AUPR);

% set test parameters
par.thresh = 0;
par.alpha = P(nP);
par.pcc = pcc;              % is the same for grnnet_eval and grnnet

% apply inference method
gsidx = correct_grn_network(grnnet,par);

% truncate predictions to 100000 interactions
gsidx = truncate_predictions(gsidx);

% evaluate prediction using various metrics
[tpr fpr prec rec l auroc aupr p_auroc p_aupr] = DREAM5_Challenge4_Evaluation(gold_edges, gsidx, pdf_aupr, pdf_auroc);

perf_total = -(log10(p_auroc)+log10(p_aupr))/2;
fprintf('\nInference method %g using only TF-NTF interactions with correction: AUPR %g   AUROC %g   Total %g\n\n',midx,aupr*100,auroc*100,perf_total)
