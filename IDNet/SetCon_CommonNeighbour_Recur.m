function [X_out,idc,p_out,recip] = SetCon_CommonNeighbour_Recur(Nsyn, Nin, Nout, pCon, pSelfCon, pRec, set_flag, X)
% Generates a random connectivity matrix X[Nout, Nin] with
% connection prob. pCon and implements the "common neighbour rule"
% (Perin et al. 2011, PNAS). 
%
% INPUT:
%   Nsyn:           Number of synapses that are already created
%   Nin:            Number of input neurons
%   Nout:           Number of output neurons
%   pCon:           Connection probability
%   pSelfCon        Connection probability for self-connections (autapses)
%   pRec            Fraction of reciprocal (bidirectional) connections
%   set_flag:       Change connection matrix according to the common
%                   neighbour rule if set_flag==TRUE
%   X:              Connection matrix between the input and output neurons
% 
% OUTPUT:
%   X_out:          Connection matrix containing idc at non-zero entries
%   idc:            Synapse indices (starting with Nsyn)
%   p_res:          Connection probability as a function of the number of
%                   common neighbours
%   recip           Actual fraction of reciprocal connections



% Set probabilities according to the number of neighbours 
% (derived from Perin et al. 2011) 
slope = 20*3.9991/Nout;
max_neigh = min(Nout,floor(1/pCon/slope));
if max_neigh==0
    disp('Warning: max_neigh=0 => CNR not applicable!');
end

% Set random connectivity with probability pCon (if not preset)
if ~exist('X', 'var')
    k1 = randperm(Nin*Nout);
    k = k1(1:round(pCon*Nin*Nout));
    X=zeros(Nin,Nout);
    idc=Nsyn+(1:length(k));
    X(k)=idc;
else
    idc = min(X(X>0)) + (1:max(X(X>0))) - 1;
    k = find(X>0)';
end;

% Compute the number of common neighbours for each pair in X
[pair_id, N_neigh] = find_neigh(X);


% Set new matrix, if desired
if ~exist('set_flag', 'var') || set_flag
        
    % Normalize pCon such that p values to add up to original pCon
    p0 = p_calc(X, pair_id, N_neigh, pCon, slope);
    p(:,1:2) = [p0(:,1), p0(:,5)];
    
    % Get bidirectional pairs
    [pairs, N_neigh_pairs] = find_neigh_recur(X);
     
    % Select connections randomly according to p (except self-connections)
    pair_id_selected = [];
    for i=1:length(p(:,1))
        pair_id_act = pair_id(N_neigh == p(i,1),:);
        pairs_act = pairs(N_neigh_pairs == p(i,1),:);
        pair_id_old = k(ismember(k, pair_id_act));
        N_rec = floor(p(i,2)*(pRec/2)*length(pair_id_act));                         % Number of recurrently connected pairs
        N_uni = ceil(p(i,2)*length(pair_id_act)) - 2*N_rec;                         % Number of unidirectionally connected pairs
        k1 = randperm(length(pairs_act(:,1)));                                      % Randomize pair indices
        
        % Determine number of unidirectional and recurrent connections in old pairs
        dummy = zeros(length(pair_id_old),2);
        dummy(:,1) = ceil(pair_id_old/Nout)';
        dummy(:,2) = pair_id_old'-(dummy(:,1)-1)*Nout;
        pair_id_old_2 = (dummy(:,2)'-1)*Nout + dummy(:,1)'; % counterpart of each pair_id in pair_id_old
        pair_id_old_uni = pair_id_old(~ismember(pair_id_old, pair_id_old_2))';
        pairs_old_rec = [pair_id_old(ismember(pair_id_old, pair_id_old_2))' pair_id_old_2(ismember(pair_id_old, pair_id_old_2))'];
        
        
        % Set bidirectional connections
        if ~isempty(pairs_old_rec) && length(pairs_old_rec(:,1)) > N_rec            % if pairs_old_rec is too long, randomly select from it
            k1_old = randperm(length(pairs_old_rec(:,1)));
            pair_id_selected = [pair_id_selected; pairs_old_rec(k1_old(1:N_rec),1)];
            pair_id_selected = [pair_id_selected; pairs_old_rec(k1_old(1:N_rec),2)];

        else                                                                        % if pairs_old_rec is too short, use random pairs from pairs_act
            pair_id_selected = [pair_id_selected; pairs_act(k1(1:N_rec),1)];        % Set bidirectional connections
            pair_id_selected = [pair_id_selected; pairs_act(k1(1:N_rec),2)];        % Set bidirectional connections (counterparts)
        end
        pair_id_act = pair_id_act(~ismember(pair_id_act,pair_id_selected));         % Do not reuse pairs that are already selected
        k1 = randperm(length(pair_id_act));                                         % Randomized pair indices
        
        % Set unidirectional connections
        if length(pair_id_old_uni) > N_uni        % if pair_id_old_uni is too long, randomly select from it
            k1_old = randperm(length(pair_id_old_uni(:,1)));
            pair_id_selected = [pair_id_selected; pair_id_old_uni(k1_old(1:N_uni))];
            
        else                                      % if pair_old is too short, use random pairs from pairs_act
            pair_id_selected = [pair_id_selected; pair_id_act(k1(1:N_uni),:)];
        end;
%        i
    end;
    
    
    % Select self-connections randomly according to pSelfCon 
    pair_id_diag = ((1:Nout)-1)*Nout + (1:Nout);
    pair_id_diag_old = find(X(pair_id_diag)>0);
    k1 = randperm(length(pair_id_diag));
    if length(pair_id_diag_old)>round(pSelfCon*length(pair_id_diag))        % if ind_diag_old is too long, randomly select from it
        pair_id_diag_old = k(ismember(k, pair_id_diag));
        k1_old = randperm(length(pair_id_diag_old));
        pair_id_selected = [pair_id_selected; pair_id_diag_old(k1_old(1:round(pSelfCon*length(pair_id_diag))))'];
        
    else                                                                    % if ind_diag_old is too short, use random pairs from pair_id_diag
        N_diag = round(pSelfCon*length(pair_id_diag))-length(pair_id_diag_old);
        pair_id_selected = [pair_id_selected; pair_id_diag(k1(1:N_diag))'];
    end;

    
    % Set connectivity matrix
    rand_ind = randperm(length(pair_id_selected));      % randomize indices
    X_out=zeros(Nin,Nout);
    idc=Nsyn+(1:length(pair_id_selected));
    X_out(pair_id_selected)=idc;
    X_out(pair_id_selected)=idc(rand_ind);
        
else
    X_out = X;
    disp('Warning: No CNR was applied! Correct probabilities for self-connections might be unregarded!');
end;


% Analyse connectivity
[pair_id_out, N_neigh_out] = find_neigh(X_out);
p_out = p_calc(X_out, pair_id_out, N_neigh_out, pCon, slope);

% Compute fraction of reciprocal connections
test = (X_out>0) + (X_out>0)';
recip = length(test(test>1)) / length(X_out(X_out>0));



% -------------------------------------------------------------------------
% ---------------------  Auxillary functions  -----------------------------
% -------------------------------------------------------------------------

function [pair_id, N_neigh] = find_neigh(X)
% function to compute the number of common neighbours for each pair of
% neurons connected by X

Nin = length(X(1,:));
Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

% Find pairs which share at least one common neighbour
X2 = (X + X')>0;             % make matrix symmetrical, as only pairs are important, not directionality of connections
dummy = cumsum(X2>0,2); 
pairs = zeros(1,2);
kk=0;
for i=2:max(max(dummy))    % loop up to maximal number of neighbours of any of the cells
    idx = find(dummy(:,end) == i);    % neurons with i outgoing connections
    neighbours = zeros(length(idx),i);
    for j=1:length(idx)
        neighbours(j,:) = find(X2(idx(j),:)>0);
        pairs(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = nchoosek(neighbours(j,:),2);     % Get all combinations of neurons that project to that neuron as pairs
        kk = kk+nchoosek(length(neighbours(j,:)),2);
    end;
end;

% Make indices for bidirectional connections for each pair
% and compute number of neighbours
pair_id = [(pairs(:,1)-1)*Nout + pairs(:,2); (pairs(:,2)-1)*Nout + pairs(:,1)];
[pair_id, ~, pair_ind_redund] = unique(sort(pair_id));
N_neigh = accumarray(pair_ind_redund,1);

% Merge with random pairs, keep only unique pairs (as neighboured)
out_ind = ~ismember(1:Nin*Nout, pair_id);
N_neigh = [zeros(length(find(out_ind)),1); N_neigh];
[pair_id, sort_idx] = sort([find(out_ind)'; pair_id]);
N_neigh = N_neigh(sort_idx);

% Set self-connections to "-1 neighbours"
N_neigh(ind_diag) = -1;



function [pair_id, N_neigh] = find_neigh_recur(X)
% function to compute the number of common neighbours for each pair of
% neurons connected by X

Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

% Find pairs which share at least one common neighbour
X2 = (X + X')>0;                            % make matrix symmetrical, as only pairs are important, not directionality of connections
dummy = cumsum(X2>0,2); 
pairs = zeros(1,2);
kk=0;
for i=2:max(max(dummy))    
    idx = find(dummy(:,end) == i);          % neurons with i outgoing connections
    neighbours = zeros(length(idx),i);
    for j=1:length(idx)
        neighbours(j,:) = find(X2(idx(j),:)>0);
        pairs(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = nchoosek(neighbours(j,:),2);
        kk = kk+nchoosek(length(neighbours(j,:)),2);
    end;
end;

% Make indices for bidirectional connections for each pair
% and compute number of neighbours
pair_id(:,1) = (pairs(:,1)-1)*Nout + pairs(:,2);             
pair_id(:,2) = (pairs(:,2)-1)*Nout + pairs(:,1);
[~, pair_ind_unique, pair_ind_redund] = unique(pair_id(:,1));
pair_id = pair_id(pair_ind_unique,:);
N_neigh = accumarray(pair_ind_redund,1);

% Get connections on the lower triangluar part of X
X_tril = tril(ones(size(X)));
tril_ind = find((X_tril>0));
out_ind = ~ismember(tril_ind, pair_id(:,1));
pair_id_zero(:,1) = tril_ind(out_ind);

% Reconstruct connections on the upper triangluar part of X
pairs_zero(:,1) = ceil(pair_id_zero(:,1)/Nout);
pairs_zero(:,2) = pair_id_zero(:,1)-(pairs_zero(:,1)-1)*Nout;
pair_id_zero(:,2) = (pairs_zero(:,2)-1)*Nout + pairs_zero(:,1);

% Add zero neigbour pairs (all the rest of X_tril)
N_neigh = [zeros(length(find(out_ind)),1); N_neigh];
pair_id = [pair_id_zero; pair_id];
[~, sort_idx] = sort(pair_id(:,1));
pair_id = pair_id(sort_idx,:);
N_neigh = N_neigh(sort_idx);

% Set self-connections to "-1 neighbours"
N_neigh(ismember(pair_id(:,1),ind_diag)) = -1;



function p = p_calc(X, pair_id, N_neigh, pCon, slope)
% function to analyse connectivity

Nout = length(X(:,1));

% Take out self-connections
ind_diag = ((1:Nout)-1)*Nout + (1:Nout);
X(ind_diag) = zeros(1,Nout);

p(:,1) = unique(N_neigh(N_neigh>=0));    % do not use self-connections
pair_id_selected = find(X>0);
N_neigh_selected = N_neigh(ismember(pair_id, pair_id_selected));
for i=1:length(p(:,1))
    pair_act = pair_id(N_neigh == p(i,1));
    p(i,2) = length(pair_act);
    p(i,3) = sum(N_neigh_selected==p(i,1));
    p(i,4) = sum(N_neigh_selected==p(i,1)) / length(pair_act);
end;

% Normalize pCon such that p values to add up to original pCon
off = fminsearch(@(N0) p_min(N0, p, pCon, slope), max(p(:,1)));
N1 = floor(min(max(p(:,1)),  off+1/slope/pCon));
N0 = max(min(p(:,1)), ceil(off));
ind = ismember(p(:,1), N0:N1);
p(:,5) = zeros(size(p(:,1)));
p(ind,5) = pCon*slope*(p(ind,1)-off);
p(find(ind,1,'last')+1:end,5) = 1;
if ~isempty(find(p(:,5)<0, 1))
    disp('Warning: negative probabilities!');
    n_idx=find(p(:,5)<0);
    disp(['i=' num2str(n_idx)]);
elseif ~isempty(find(p(:,5)>1, 1))
    disp('Warning: probabilities > 1');
    p_idx=find(p(:,5)<0);
    disp(['i=' num2str(p_idx)]);
end


function res = p_min(N0, p, pCon, slope)
% function to be minimized to compute offset for p_Neighbours

N1 = round(min(max(p(:,1)),  N0+1/slope/pCon));

N_2 = p(p(:,1)>N0 & p(:,1)<=N1,1);
pN_2 = p(p(:,1)>N0 & p(:,1)<=N1,2)/sum(p(:,2));
pN_3 = p(p(:,1)>N1,2)/sum(p(:,2));

res = abs(sum([pN_2*slope.*(N_2-N0); pN_3]) -1);


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim