function [idx,X] = get_idx(M,nSites)
% GET_IDX  Concise index sets for a modularization.
%    IDX = GET_IDX(M,NSITES) for modularization M and total number of observable components NSITES
%    (without the reference variable s_ref) returns the indices of the elements of the 
%    NSITES-by-NSITES moment ratio matrix (obtained by GET_B) that correspond to concise index sets.
%
%    Each column of IDX contains the indices of a concise index set, where 0 denotes no index.
%
%    IDX is an input argument for the function GET_T.
%
% From: "Correlations reveal the hierarchical organization of networks with latent binary variables" (2023) Stefan Häusler

% (c) 2023 Stefan Häusler
% This code is licensed under BSD-3-Clause license (see LICENSE for details)



% check arguments

if ~iscell(M)
   error('M must be of type CELL.')
end

if ~isnumeric(nSites)|(prod(size(nSites))~=1)|(nSites<1)
   error('nSites must be a positive scalar.')
end


if isempty(M)
   error('Cells of M must contain nonempty vectors.')
end

nM = length(M);

for i = 1:nM
  if isempty(M{i})|(length(M{i})~=prod(size(M{i})))
    error('Cells of M must contain nonempty vectors.')
  end
  if any(M{i}>nSites)
    error('nSites too small.')
  end
  if length(M{i})==1
    error('Cells of M must not contain scalars.')
  end
end

% extend modularization

for i = 1:nSites
  M{end+1} = i;
end
nM = length(M);

% generate redundant index sets

Y = [];
for i = 1:nM
  for j = (i+1):nM
    if isempty(intersect(M{i},M{j}))
      b = zeros(nSites);
      b(M{i},M{j}) = 1;
      b(M{j},M{i}) = 1;
      b = triu(b,1);
      Y{end+1}=find(b)';
  
    elseif ~( all(ismember(M{i},M{j}))|all(ismember(M{j},M{i})) )
      error('No hierarchical modularization.')
    end  
  end
end

% generate concise index sets

nY = length(Y);
for i = nY:-1:1
  c = 0; 
  for j = 1:i-1
%    if all(ismember(Y{i},Y{j})) % too slow
    if (length(Y{i})<length(Y{j})) && any(Y{i}(1)==Y{j})  % faster version of line above (works due to partitions)
      c = 1;
    end
  end
  if c
    Y(i) = [];
  end  
end
X = Y;
nX = length(X);


% generate idx for get_T

idx = [];
for i = 1:nX
  if length(X{i})>1
    idx(1:length(X{i}),end+1) = X{i};
  end
end


