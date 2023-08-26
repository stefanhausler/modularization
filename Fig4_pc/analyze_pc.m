% Statistical test for functional modularizations of a CA1 pyramidal cell model.
%
% Reproduces results shown in Fig 4 discussed in the section "Inferring flat modularizations in dendrites of pyramidal neurons" in
% "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% Runtime on an Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (28 cores):  < 60 sec execution time and < 3 GB memory
%
% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)


clear all

addpath('../lib')    % path to test library


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% data parameters

Tstim = 20;          % stimulation duration 20 or 60 [min]
datadir = 'data';    % data directory

% test parameters

par.thresh = 5;      % threshold of 5, 6 or 7
par.METHOD = 2;      % 0 ... uncorrected, 1 ... simple correction, 2 ... advanced correction

% all single functional modules (module -by- network component)

S = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];

% all single functional modules to be combined and tested for linearity (module -by- network component)

SL = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
      0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];

% generate all combinations of modules in SL with 11 sites

nSL = size(SL,1);
nSites = size(SL,2);
SC = zeros(2^nSL-1,nSites);             % number of combinations without restriction to 11 sites
for j = 1:size(SC,1)
    it = (dec2bin(j,nSL))==49;
    SC(j,:) = (sum(SL(it,:),1) > 0);
end
SC = unique(SC,'rows');
SC(sum(SC,2)~=11,:) = [];               % remove all combinations with more or less than 11 sites

% additional modularizations to be tested

M1 = {[5 6],[12:18],[1:4 7:11 19:26]};  % modules 4, 7 and remaining sites (nr according to the paper)
M2 = {[4:7],[12:18],[1:3 8:11 19:26]};  % modules 2, 7 and remaining sites (nr according to the paper)
S0 = 1:nSites;                          % M_L

% initialize concice index sets of all modularizations

idx0 = get_idx_linear(S0,nSites); % indices for concise index sets (columns)
idx1 = get_idx(M1,nSites);              % indices for concise index sets (columns)
idx2 = get_idx(M2,nSites);              % indices for concise index sets (columns)

for j = 1:size(S,1)
    fprintf('Initialize concise index set for single module %g/%g                  \r',j,size(S,1))
    M = {find(S(j,:)==1)};               % modularization is j-th single functional module
    IDX{j} = get_idx(M,nSites);          % indices for concise index sets (columns)
end % j
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statistical test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for nTr = 1:10;
    
    fname = sprintf('%s/data_recording_PC_%gmin_trial%g.mat',datadir,Tstim,nTr);
    fprintf('Load %s\n',fname)
    load(fname,'Vm','VmSites')     % load data file
    % Vm      ... 1 -by- time point
    % VmSites ... recording site -by- time point
    
    nSites = size(VmSites,1);      % number of recording sites
    nS = size(VmSites,2);          % number of time points
    
    stat = get_B(Vm,VmSites,par);  % get moment ratios, variances, covariance matrix and indices of elememts to clamp
    
    b(:,:,nTr) = stat.b;           % store for visualization
    
    for j = 1:size(S,1)
        fprintf('Single module %g/%g                  \r',j,size(S,1))
        
        p(j,nTr) = get_T(stat,IDX{j},par);    % calculate statistic and its p-value
        
    end % j
    
    for j = 1:size(SC,1)
        fprintf('Combination of modules %g/%g                  \r',j,size(SC,1))
        
        S1 = find(SC(j,:));                 % single function module containing 11 components tested for linearity
        idx = get_idx_linear(S1,nSites);    % indices for concise index sets (columns)
        pL(j,nTr) = get_T(stat,idx,par);    % calculate statistic and its p-value
        
    end % j
    
    % test additional modularizations
    p0 = get_T(stat,idx0,par);              % calculate statistic and its p-value
    p1 = get_T(stat,idx1,par);              % calculate statistic and its p-value
    p2 = get_T(stat,idx2,par);              % calculate statistic and its p-value
    
    P(:,nTr) = [log10([p0 p1 p2])];         % store results
    
end % nTr
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% determine the maximum p-value of all combinations in SL that contain the single functional module [1 2 5:size(SL,1)] in SL
clear x
for n = [1 2 5:size(SL,1)]
    i = find(SC*SL(n,:)' == sum(SL(n,:)));
    x(n,:) = max(pL(i,:),[],1);
end

% determine the maximum p-value of all combinations in SL that contain the single functional module [3 4] in SL
for n = [3 4]
    i = find( (SC*SL(n,:)' == sum(SL(n,:))) & ~(SC*SL(2,:)' == sum(SL(2,:))));
    x(n,:) = max(pL(i,:),[],1);
end

% calculate significance levels

amin = get_amin(par,nSites); % minimum significance level
alpha = 0.05;  % nominal significance level

nH = 15+2+1; % number of hypotheses = 15 + (modules 2&7) + (modules 4&7) + linear
a = (alpha-amin)/(1-amin);
alpha1 = a/(2*nH); % /2 correction, because of multiple testing in two panels (Bonferroni correction)

nH = size(pL,1); % number of hypotheses = all combinations in SL
a = (alpha-amin)/(1-amin);
alpha2 = a/(2*nH); % /2 correction, because of multiple testing in two panels (Bonferroni correction)

disp('Significance of three extra modularizations M_L, M1 and M2')
disp(median(P'))
disp(median(P')<log10(alpha1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure

subplot(3,2,5,'align');

B = mean(b,3);
B = B + diag(nan(1,length(B)));
caxis([-3.9 3.3]);
imagesc(B)
colormap(jet)
title('B')
cbh = colorbar;
set(cbh,'ytick',[-3.9 3.3])
set(cbh,'yticklabel',[-3.9 3.3])
ylim(cbh,[-3.9 3.3])


subplot(4,3,[1 2],'align');
hold on

if Tstim == 60
    marker = 'o';
    coloridx = 5;
    color = 0.5*[1 1 1];
    yl = 10^-120;
else
    marker = 'x';
    coloridx = 4;
    color = [0 0 0];
    yl = 10^-30;
end

plot(squeeze(p),'color',color,'marker',marker,'markersize',3,'linestyle','none');
plot([[1:15]-0.35; [1:15]+0.35],repmat(10.^squeeze(median(log10(p),2))',[2 1]),'color',color,'marker','none');

set(gca,'yscale','log')
box off
axis tight
ylim([1e-20 1])
set(gca,'ytick',[1e-20 1e-10 1])
xlim([0 16])
set(gca,'XAxisLocation','top')
hold on
plot(xlim,[1 1]*alpha1,'k:')
set(gca,'xtick',[1 15])
set(gca,'xticklabel',{'S_1' 'S_{15}'})
ylabel('p-value')

subplot(4,3,[4 5],'align');
hold on

i = [1 2 4:13]; % to simplify panel don't show results for S4
xc = [1:2:23]';

plot(repmat(xc,[1,10]),x(i,:),'markersize',3,'color',color,'marker',marker,'linestyle','none')
plot([xc-0.35, xc+0.35]',repmat(10.^median(log10(x(i,:)),2)',[2 1]),'color',color,'marker','none');

set(gca,'xtick',xc)
set(gca,'xticklabel',{'S_1' 'S_2' 'C_{ }' 'S_5' 'S_7' 'S_{13}' 'S_{14}' 'S_{15}' 's_{3}' 's_{11}' 's_{19}' 's_{22}' })
set(gca,'XAxisLocation','top')
set(gca,'yscale','log')
ylim([yl 1])
grid on
set(gca,'ytick',[1e-120 1e-90 1e-60 1e-30 1])
xlim([0 2*size(SL,1)-1])
plot(xlim,[1 1]*alpha2,'k:')
grid off
box off
ylabel('p-value')

drawnow
