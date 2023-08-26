% Statistical test for functional modularizations of a neural network.
%
% Reproduces the results shown in Fig. 5 discussed in the section "Inferring nested modularizations in neural networks" in
% "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% Runtime on an Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (28 cores):  < 40 sec execution time and < 1 GB memory
%
% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)


clear all

addpath('../lib')    % path to test library


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% data parameters

datadir = 'data';                     % data directory
datafile = 'data_recording_nn.mat';   % data filename
% datafile = 'data_recording_nn_linear.mat';   % second data filename

% test parameters

par.thresh = 5;      % threshold of 5, 6 or 7
par.METHOD = 2;      % 0 ... uncorrected, 1 ... simple correction, 2 ... advanced correction

% load data

fname = sprintf('%s/%s',datadir,datafile);
fprintf('Load %s\n',fname)
load(fname,'x','z')     
% x ... 4 -by- no. of time intervals -by- no. of trials 
% z ... 1 -by- no. of time intervals -by- no. of trial 

nSites = size(x,1);     % number of network components (without z)
nS  = size(x,2);        % number of time interval
nTr = size(x,3);        % number of trials (in the article 1e5)

disp('Normalize to mean zero and unit variance.')

x = x - repmat(mean(x,2),[1 nS 1]);
x = x./repmat(std(x,[],2),[1 nS 1]);
z = z - repmat(mean(z,2),[1 nS 1]);
z = z./repmat(std(z,[],2),[1 nS 1]);

disp('Calculate mean correlation matrix averaged over trials.')

CC = zeros(5);
for i = 1:nTr
    CC = CC + corrcoef(cat(1,x(:,:,i),z(:,:,i))');
end
CC = CC/nTr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statistical test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Generate modularizations.')

% All 25 modularizations (excl. M_L), see Figure S3 in the paper
M = {{[3,4]},
    {[2,4]},
    {[2,3]},
    {[1,4]},
    {[1,3]},
    {[1,2]},
    {[1,2],[3,4]},
    {[1,3],[2,4]},
    {[1,4],[2,3]},
    {[2,3,4]},
    {[1,3,4]},
    {[1,2,4]},
    {[1,2,3]},
    {[2,3,4],[3,4]},
    {[3,2,4],[2,4]},
    {[4,2,3],[2,3]},
    {[1,3,4],[3,4]},
    {[3,1,4],[1,4]},
    {[4,1,3],[1,3]},
    {[1,2,4],[2,4]},
    {[2,1,4],[1,4]},
    {[4,1,2],[1,2]},
    {[1,2,3],[2,3]},
    {[2,1,3],[1,3]},
    {[3,1,2],[1,2]}};

nM = length(M);     % number of modularizations (excl. M_L)

disp('Generate concise index sets.')

IDX = [];
for i = 1:nM
    fprintf('Initialize concise index set for modularization %g/%g                  \r',i,nM)
    IDX{i} = get_idx(M{i},nSites);          % indices for concise index sets for all modularizations except M_L
end % i
fprintf('\n')

IDXL = get_idx_linear(1:nSites,nSites);     % indices for concise index sets of M_L

% Evaluate modularizations

p = [];
for j = 1:nTr                               % loop over trials
    fprintf('Evaluate modularizations (%g%%)                \r',j/nTr*100)
    stat = get_B(z(:,:,j),x(:,:,j),par);    % get moment ratios, variances, covariance matrix and indices of elememts to clamp
    for i = 1:nM                            % loop over modularizations 
        p(i,j) = get_T(stat,IDX{i},par);    % calculate statistic and its p-value
    end % i
    p(nM+1,j) = get_T(stat,IDXL,par);       % calculate statistic and its p-value for M_L
end % j
fprintf('\n')

logit_p = log10(p)';                        % transform to log scale

disp('Calculate Bonferroni corrected individual significance level.')

amin = get_amin(par,nSites); % minimum significance level
alpha = 0.01;                % nominal significance level
nH = nM+1;                   % number of hypotheses
a = (alpha-amin)/(1-amin);   % individual significance level
alpha_corr = a/nH;           % Bonferroni corrected individual significance level


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
clf

% plot correlation matrix

subplot(4,4,1,'align')
cla

imagesc(CC-diag(diag(CC)))

p = get(gca,'position');
cbh = colorbar;
set(cbh,'ytick',[0 0.9])
ylabel(cbh,'cc')
ylim(cbh,[0 0.9])
caxis([0 0.9])

set(gca,'position',p)
set(gca,'xtick',[1:5])
set(gca,'xticklabel',{'s_1','s_2','s_3','s_4','s_5'})
set(gca,'ytick',[1:5])
set(gca,'yticklabel',{'s_1','s_2','s_3','s_4','s_5'})
set(gca,'TickLength',[0 0])

% plot p-values

subplot(3,1,2,'align')
cla
hold on

% Calculate 95 perc. confidence interval

l_p = prctile(logit_p,0.025*100);
m_p = mean(logit_p);
h_p = prctile(logit_p,(1-0.025)*100);

% order of modularizations 

k = [1 10 14 15 16 2 3 20 23 5 4 17 18 19 12 13 6 7 9 24 25 8 21 22 11 26];

% plot results

errorbar(1:26,m_p(k),l_p(k)-m_p(k),h_p(k)-m_p(k),'marker','none','linestyle','none','color',[0 0 0])
boxplot(logit_p(:,k),'Notch','off','color',[0 0 0],'symbol','','whisker',0);
plot(xlim,log10(alpha_corr)*[1 1],'k--')

xlim([0.5 26.5])
ylim([-32 0])

set(gca,'ytick',[-32 -16 0])
set(gca,'yticklabel',{1e-32 1e-16 1e-0'})

set(gca,'xtick',[1:26])
set(gca,'xticklabel',[k 0])
set(gca,'XAxisLocation','top');

box off
xlabel('Modularization index')
ylabel('p-value')

drawnow



