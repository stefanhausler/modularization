function [x,info] = get_grn_data(NETWORK,path_dream5)
% GET_GRN_DATA  Load gene regulatory network of the DREAM5 challenge.
%    [X,INFO] = GET_GRN_DATA(NETWORK,PATH_DREAM5) for path PATH_DREAM5 and index NETWORK returns
%    gene expression levels X (#genes-by-#samples) and network information
%
%       INFO.tfidx   indices of all transcription factors (#(transcription factors)-by-1)
%       INFO.gsidx   gold standard (#interactions-by-3)
%
%    PATH_DREAM5 (string) specifies the location of the DREAM5 folder.
%
%    NETWORK index 1 ... In silico network
%                  3 ... E. coli network
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)


info = [];

st = dbstack;
datapath = fileparts(which(st(1).name));

if exist(sprintf('%s/dataset%d.mat',datapath,NETWORK))
    load(sprintf('%s/dataset%d.mat',datapath,NETWORK))
else
    
    % Network1: In silico / 4012 interactions
    % Network2: S. aureus (not used for scoring) / 5018 interactions
    % Network3: E. coli / 2066 interactions
    % Network4: S. cerevisiae / 3940 interactions
    datasets = {'Network1','Network2','Network3','Network4'};
    
    tb = readtable( sprintf('%s/%s/input_data/net%d_expression_data.csv',path_dream5,datasets{NETWORK},NETWORK),'delimiter','tab');
    tf = readtable( sprintf('%s/%s/input_data/net%d_transcription_factors.csv',path_dream5,datasets{NETWORK},NETWORK),'delimiter','tab','ReadVariableNames',0);
    gs = readtable( sprintf('%s/%s/gold_standard/DREAM5_NetworkInference_GoldStandard_Network%d.csv',path_dream5,datasets{NETWORK},NETWORK),'delimiter','tab','ReadVariableNames',0);
    
    tfidx = zeros(size(tf,1),3);;
    for i = 1:size(tf,1)
        tfidx(i) = str2num(tf{i,1}{1}(2:end));
    end
    
    gsidx = zeros(size(gs,1),3);
    for i = 1:size(gs,1)
        fprintf('%g               \r', i/size(gs,1))
        gsidx(i,1) = str2num(gs{i,1}{1}(2:end));
        gsidx(i,2) = str2num(gs{i,2}{1}(2:end));
        gsidx(i,3) = gs{i,3};
    end
    
    x = tb{1:end,1:end}';
    VariableNames = tb.Properties.VariableNames;
    
    vn = [];
    for i = 1:length(VariableNames)
        vn(i) = str2num(VariableNames{i}(2:end));
    end
    
    if ~all(vn==[1:length(vn)])
        error('Check indices!')
    end
    
    for i = 1:size(x)
        x(i,:) = x(i,:)-median(x(i,:));
        x(i,:) = x(i,:)/std(x(i,:));
    end
    
    info.tfidx = tfidx;
    info.gsidx = gsidx;
    
    save(sprintf('%s/dataset%d.mat',datapath,NETWORK),'x','info')
end

end

