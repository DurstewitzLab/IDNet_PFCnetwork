function [X_out,idc,p_res] = SetCon_CommonNeighbour_Cross(Nsyn, X_11,X_22,pCon,para_p_input,para_p_output,set_flag,X)
% Generates a random connectivity matrix X[Nout, Nin] across two different
% neuron populations with connection probabilities pCon and implements the 
% "common neighbour rule" (Perin et al. 2011, PNAS). 
%
% INPUT:
%   Nsyn:           Number of synapses that are already created
%   X11:            Connection matrix within the output neurons
%   X22:            Connection matrix within the input neurons
%   pCon:           Connection probability
%   para_p_input:   Common neighbour rule parameters for the input neurons
%   para_p_output:  Common neighbour rule parameters for the output neurons
%   set_flag:       Change connection matrix according to the common
%                   neighbour rule if set_flag==TRUE
%   X:              Connection matrix between the input and output neurons
% 
% OUTPUT:
%   X_out:          Connection matrix containing idc at non-zero entries
%   idc:            Synapse indices (starting with Nsyn)
%   p_res:          Connection probability as a function of the number of
%                   common neighbours


Nin = length(X_22(1,:));   
Nout = length(X_11(:,1));

% Set probabilities according to the number of neighbours 
% (derived from Kampa, Letzkus and Stuart, 2006) 
if ~exist('para_p_input', 'var') || isempty(para_p_input')
    para_p_input = [0.1411   0.1677];
end
if ~exist('para_p_output', 'var')|| isempty(para_p_output')
    para_p_output = [-0.1621   0.3514];
end
max_neigh_input = floor((1/pCon-para_p_input(2))/para_p_input(1));              % ensures p <= 1
p_Neighbours_input = para_p_input(1)*(1:max_neigh_input) + para_p_input(2);
max_neigh_output = floor(-para_p_output(2)/para_p_output(1));                   % ensures p >= 0
p_Neighbours_output = para_p_output(1)*(1:max_neigh_output) + para_p_output(2);

% Set random connectivity with probability pCon (if not preset)
if ~exist('X', 'var')
    k1 = randperm(Nin*Nout);
    k = k1(1:round(pCon*Nin*Nout));
    X=zeros(Nin,Nout);
    idc=Nsyn+(1:length(k));
    X(k)=idc;
else
    idc = min(X(X>0)) + (1:max(X(X>0))) - 1;
    k = X(X>0)';
end;

% Compute the number of common neighbours for each pair in X
[~,       ~, neigh_00] = find_neigh(X,X_22,0,0);  % neurons projecting to a common neighbour in the output layer
[~,       ~, neigh_01] = find_neigh(X,X_22,0,1);  % the same, but the  output neuron in the pair receives input from the neighbour instead of projecting to it
[~,       ~, neigh_10] = find_neigh(X,X_11,1,0);  % neurons receiving input from a common neighbour in the input layer
[pair_id, ~, neigh_11] = find_neigh(X,X_11,1,1);  % the same, but the input neuron in the pair projects to the neighbour instead of receiving input from it

N_neigh_input = zeros(length(pair_id),1);
N_neigh_output = zeros(length(pair_id),1);
for i=1:length(pair_id)
    neigh_act = unique([neigh_00(i,:), neigh_01(i,:)]);
    N_neigh_output(i) = length(neigh_act(neigh_act>0));
    neigh_act = unique([neigh_10(i,:), neigh_11(i,:)]);
    N_neigh_input(i) = length(neigh_act(neigh_act>0));
end;


% Set new matrix, if desired
if ~exist('set_flag', 'var') || set_flag
        
    % Normalize pCon such that p values to add up to original pCon
    p = p_calc_cross(X, pair_id, N_neigh_input, N_neigh_output, pCon, p_Neighbours_output'*p_Neighbours_input);
    
    % Select connections randomly according to p
    pair_id_selected = [];
    for i=1:length(p(:,1,1))
        for j=1:length(p(1,:,2))
            pair_act = pair_id(N_neigh_output == p(i,j,1) & N_neigh_input == p(i,j,2));
            pair_old = k(ismember(k, pair_act));
            
            if length(pair_old) > round(p(i,j,6)*length(pair_act))        % if pair_old is too long, randomly select from it
                k1 = randperm(length(pair_old));
                pair_id_selected = sort([pair_id_selected, pair_old(k1(1:round(p(i,j,6)*length(pair_act))))]);                  % apply common neighbour rule
                
            else                                                        % if pair_old is too short, use all pairs in it and add random pairs from pair_act
                k1 = randperm(length(pair_act));
                pair_id_selected = sort([pair_id_selected, pair_act(k1(1:round(p(i,j,6)*length(pair_act))))']);    % apply common neighbour rule
            end;
        end;
    end;

    % Set connectivity matrix
    rand_ind = randperm(length(pair_id_selected));      % randomize indices
    X_out=zeros(Nout,Nin);
    idc=Nsyn+(1:length(pair_id_selected));
    X_out(pair_id_selected)=idc;
    X_out(pair_id_selected)=idc(rand_ind);
        
else
    X_out = X;
end;


% Analyse connectivity
[~,           ~, neigh_00_res] = find_neigh(X,X_22,0,0);  % neurons projecting to a common neighbour in the output layer
[~,           ~, neigh_01_res] = find_neigh(X,X_22,0,1);  % the same, but the  output neuron in the pair receives input from the neighbour instead of projecting to it
[~,           ~, neigh_10_res] = find_neigh(X,X_11,1,0);  % neurons receiving input from a common neighbour in the input layer
[pair_id_res, ~, neigh_11_res] = find_neigh(X,X_11,1,1);  % the same, but the input neuron in the pair projects to the neighbour instead of receiving input from it

N_neigh_input_res = zeros(length(pair_id_res),1);
N_neigh_output_res = zeros(length(pair_id_res),1);
for i=1:length(pair_id_res)
    neigh_act = unique([neigh_00_res(i,:), neigh_01_res(i,:)]);
    N_neigh_output_res(i) = length(neigh_act(neigh_act>0));
    neigh_act = unique([neigh_10_res(i), neigh_11_res(i)]);
    N_neigh_input_res(i) = length(neigh_act(neigh_act>0));
end;
p_res = p_calc_cross(X, pair_id, N_neigh_input_res, N_neigh_output_res, pCon, p_Neighbours_output'*p_Neighbours_input);



% -------------------------------------------------------------------------
% ---------------------  Auxillary functions  -----------------------------
% -------------------------------------------------------------------------

function [pair_id, N_neigh, neigh_out] = find_neigh(X,X_rec,in_flag,switch_flag)
% function to compute the number of common neighbours for each pair of
% neurons connected by X 

% in_flag = 1: common input neurons in X_rec (output neurons otherwise)
% switch_flag = 1: direction of X connections is reversed

Nin = length(X(1,:));   
Nout = length(X(:,1));

% Find pairs which share at least one common neighbour
if in_flag == 1
    if switch_flag == 1
        X_act = [X' X_rec]>0;
        dim = 2;
    else 
        X_act = [X; X_rec]>0;
        dim = 1;
    end;
else
    if switch_flag == 1
        X_act = [X'; X_rec]>0;
        dim = 1;
    else 
        X_act = [X X_rec]>0;
        dim = 2;
    end;
end;
dummy = cumsum(X_act,dim);
pairs = zeros(1,2);
kk=0;

neigh=[];
for i=2:max(max(dummy))    % loop up to maximal number of neighbours of any of the cells
    if dim == 1
        idx = find(dummy(end,:) == i);    % neurons with i outgoing connections
    else
        idx = find(dummy(:,end) == i);    % neurons with i incoming connections
    end;
    neighbours = zeros(length(idx),i);
    for j=1:length(idx)
        if dim == 1
            neighbours(j,:) = find(X_act(:,idx(j))>0);
        else
            neighbours(j,:) = find(X_act(idx(j),:)>0);
        end;
        pairs(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = nchoosek(neighbours(j,:),2);
        neigh(kk+1:kk+nchoosek(length(neighbours(j,:)),2),:) = idx(j);
        kk = kk+nchoosek(length(neighbours(j,:)),2);
    end;
end;

% Use only those pairs which cross the layers
if in_flag == 1
    if switch_flag == 1
        neigh = neigh((pairs(:,1)<=Nout & pairs(:,2)>Nout) | (pairs(:,1)>Nout & pairs(:,2)<=Nout));
        pairs = pairs((pairs(:,1)<=Nout & pairs(:,2)>Nout) | (pairs(:,1)>Nout & pairs(:,2)<=Nout),:);
    else
        neigh = neigh((pairs(:,1)<=Nout & pairs(:,2)>Nout) | (pairs(:,1)>Nout & pairs(:,2)<=Nout));
        pairs = pairs((pairs(:,1)<=Nout & pairs(:,2)>Nout) | (pairs(:,1)>Nout & pairs(:,2)<=Nout),:);
    end;
else
    if switch_flag == 1
        neigh = neigh((pairs(:,1)<=Nin & pairs(:,2)>Nin) | (pairs(:,1)>Nin & pairs(:,2)<=Nin));
        pairs = pairs((pairs(:,1)<=Nin & pairs(:,2)>Nin) | (pairs(:,1)>Nin & pairs(:,2)<=Nin),:);
    else
        neigh = neigh((pairs(:,1)<=Nin & pairs(:,2)>Nin) | (pairs(:,1)>Nin & pairs(:,2)<=Nin));
        pairs = pairs((pairs(:,1)<=Nin & pairs(:,2)>Nin) | (pairs(:,1)>Nin & pairs(:,2)<=Nin),:);
    end;
end;

% Transform neuron indicies and make pair indicies
trafo = 1:Nin+Nout;
if in_flag == 1
    if switch_flag == 1
        trafo((1:Nin)+Nout) = 1:Nin;
        pairs(:,2) = trafo(pairs(:,2))';
        pair_id = (pairs(:,2)-1)*Nout + pairs(:,1);
    else
        trafo((1:Nin)+Nout) = 1:Nin;
        pairs(:,2) = trafo(pairs(:,2))';
        pair_id = (pairs(:,2)-1)*Nout + pairs(:,1);
    end;
else
    if switch_flag == 1
        trafo((1:Nout)+Nin) = 1:Nout;
        pairs(:,2) = trafo(pairs(:,2))';
        pair_id = (pairs(:,1)-1)*Nout + pairs(:,2);
    else
        trafo((1:Nout)+Nin) = 1:Nout;
        pairs(:,2) = trafo(pairs(:,2))';
        pair_id = (pairs(:,1)-1)*Nout + pairs(:,2);
    end;
end;

% Compute number of neighbours
[pair_id, ~, pair_ind_redund] = unique(pair_id);
N_neigh = accumarray(pair_ind_redund,1);

% Determine neighbours for each pair_id
neigh_out = zeros(length(pair_id), max(N_neigh));
for i=1:length(pair_id)
    neigh_out(i,1:N_neigh(i)) = neigh(pair_ind_redund==i);
end;

% Merge with random pairs
out_ind = ~ismember(1:(Nin*Nout), pair_id);
N_neigh = [zeros(length(find(out_ind)),1); N_neigh];
neigh_out = [zeros(length(find(out_ind)),max(N_neigh)); neigh_out];
[pair_id, sort_idx] = sort([find(out_ind)'; pair_id]);
N_neigh = N_neigh(sort_idx);
neigh_out = neigh_out(sort_idx,:);



function p = p_calc_cross(X, pair_id, N_neigh_input, N_neigh_output, pCon, p_Neighbours)
% function to analyse connectivity

pair_id_selected = find(X>0);
N_neigh_input_unique = unique(N_neigh_input);
N_neigh_output_unique = unique(N_neigh_output);
N_neigh_input_selected = N_neigh_input(ismember(pair_id, pair_id_selected));
N_neigh_output_selected = N_neigh_output(ismember(pair_id, pair_id_selected));

p = zeros(length(N_neigh_output_unique), length(N_neigh_input_unique), 6);
for i=1:length(N_neigh_output_unique)
    for j=1:length(N_neigh_input_unique)
        pair_act = pair_id(N_neigh_output == N_neigh_output_unique(i) & N_neigh_input == N_neigh_input_unique(j));
        p(i,j,1) = N_neigh_output_unique(i);
        p(i,j,2) = N_neigh_input_unique(j);
        p(i,j,3) = length(pair_act);
        p(i,j,4) = sum((N_neigh_output_selected==N_neigh_output_unique(i)) & (N_neigh_input_selected==N_neigh_input_unique(j)));
        if ~isempty(pair_act)
            p(i,j,5) = sum((N_neigh_output_selected==N_neigh_output_unique(i)) & (N_neigh_input_selected==N_neigh_input_unique(j))) / length(pair_act);
        end;
        if N_neigh_output_unique(i)+1>length(p_Neighbours(:,1)) || N_neigh_input_unique(j)+1>length(p_Neighbours(1,:))
            p(i,j,6) = 0;
        else
            p(i,j,6) = pCon*p_Neighbours(N_neigh_output_unique(i)+1,  N_neigh_input_unique(j)+1);
        end;
    end;
end;

% Normalize pCon such that p values to add up to original pCon
if isempty(p_Neighbours)
    NN = 0;
    MM = 0;
else
    NN = min(length(p_Neighbours(:,1)), length(p(:,1,3)));
    MM = min(length(p_Neighbours(1,:)), length(p(1,:,3)));
end;
p(:,:,6) = p(:,:,6) / pCon;
pCon = pCon *  pCon / sum(sum(pCon*p_Neighbours(1:NN,1:MM) .* p(1:NN,1:MM,3) / sum(sum(p(1:NN,1:MM,3))) ));
p(:,:,6) = min(p(:,:,6) * pCon);


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
