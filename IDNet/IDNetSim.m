function [STMtx,T,V,dV,Ninp]=IDNetSim(SimPar,ConMtx)
% MATLAB wrapper for C program IDNet.c to simulate abritrary biological
% neural networks. The neuron and synapse models are specified in IDNet.c,
% the parameters SimPar in ConfigIDNet.
%
% INPUT:
%   SimPar: A set of parameters (neurons, synapses, control parameters)
%   ConMtx (optionally): A network structure (random connectivity specified 
%                        in SimPar is used if this is empty).
% 
% OUTPUT:
%   STMtx: Cell array of spike times t_sp. Each cell contains t_sp of one
%          neuron
%   T:     Vector of simulation time steps for visualizing voltage traces
%   V:     Cell array of voltage traces for each neuron in the viewlist
%   Ninp:  Number of simulated neurons



% ------------------ Read parameters from SimPar --------------------------

CtrPar = SimPar.CtrPar;                             % control parameters

Nstripes = SimPar.Nstripes;                         % number of stripes
NTypes = SimPar.NTypes;                             % number of neurons of each type
InpNTypes = SimPar.InpNTypes;                       % number of input neurons of each type

CellPar = SimPar.CellPar;                           % (mean) neuron parameters
k_trans = SimPar.k_trans;                           % parameter for inv_transform_distribution
ParCov = SimPar.ParCov;                             % Cell-Array with 14 covariation matrices
N_sig = SimPar.N_sig;                               % SDs of neuron parameters
N_max = SimPar.N_max;                               % maxima of neuron parameters
N_min = SimPar.N_min;                               % minima of neuron parameters

NPList = SimPar.NPList;                             % type of each neuron

STypPar = SimPar.STypPar;                           % synapse parameters

ConPar= SimPar.ConPar;                              % (mean) connection parameters (weights etc.)
pCon = SimPar.pCon;                                 % connectivity (for random networks)
cluster_flag = SimPar.cluster_flag;                 % determines whether to use a common neighbour rule for connectivity
S_sig = SimPar.S_sig;                               % SDs of connection parameters
S_max = SimPar.S_max;                               % maxima of connection parameters
S_min = SimPar.S_min;                               % minima of connection parameters
ConParStripes = SimPar.ConParStripes;               % parameters for inter-stripe connectivity

ViewList = SimPar.ViewList;                         % neurons to be recorded in detail

InpSTtrains = SimPar.InpSTtrains;                   % input spike trains
NoisePar = SimPar.NoisePar;                         % parameters for random background activity (noise)
NoiseDistr = SimPar.NoiseDistr;                     % Distribution of random background activity
V0Par = SimPar.V0Par;                               % parameters for initial conditions

CA = SimPar.CA;                                     % cell assemblies

NeuronGroupsSaveArray = SimPar.NeuronGroupsSaveArray;   % Array specifying neuron groups to save synapse output for

UniqueNum = SimPar.UniqueNum;
% "UniqueNum" is an integer added to IDNet's temporary data files. This
% allows differentiation of files if multiple instances of IDNet are
% desired.


% -------------------- Set overall parameters -----------------------------
% Seed random number generator to a given state
% mtstream = RandStream.getGlobalStream;
% mtstream.State = SimPar.RState;
% RandStream.setGlobalStream(mtstream);

% Set network parameters
N = sum(NTypes);                % # neurons (per stripe)
NTypesN = length(NTypes);       % # neuron types
M = sum(InpNTypes);             % # input neurons (per stripe)
InpNTypesN = length(InpNTypes); % # input neuron types

% Initalize neuron and synapse parameters
NeuPar = zeros(length(CellPar(:,1)),Nstripes*N);
V0 = zeros(size(V0Par,2),Nstripes*N);
SPMtx=zeros(Nstripes*N,Nstripes*N, max(max(ConPar{1})));

NN = 0;
SynPar=[];
Nsyn=0;
dtax_back = [];

ind1 = 1:9;
ind3 = [1 2 3 4 5 6 8 9 10];


% -------------- Set neuron parameters and most synaptic connections -----------------

% loop over stripes/columns
reorder = cell(NTypesN,Nstripes);
for ii=1:Nstripes
    
    % loop over neuron types
    t_lat = zeros(NTypesN,1);
    t_lat_LIF = zeros(NTypesN,1);
    for i=1:NTypesN
        
        % Set neuron parameters and initial conditions (one individual set for each neuron)
        iset = (1:NTypes(i))+NN;
        reorder{i,ii} = iset;
        ind_out = iset;
        iii=0;
        while(~isempty(ind_out) && iii<1000)           % comment out this whole block before using update_inv_con_PSP

            % draw random numbers from a multivariate normal distribution
            NeuPar_multi=mvnrnd( zeros(1,length(ParCov{i}(:,1))),ParCov{i},length(ind_out))';

            % invert transformation of marginal distributions
            for j=ind1
                NeuPar_multi(j,:) = inv_transform_distribution2(NeuPar_multi(j,:),k_trans(j,i),CellPar(j+1,i),N_sig(j+1,i),N_min(ind3(j)+1,i));
            end
            
            % complete parameter matrix by "C" and "a"
            for j=ind1
                NeuPar(ind3(j),ind_out) = NeuPar_multi(j,:);
            end
            NeuPar(1,ind_out) = NeuPar(1,ind_out).*NeuPar(2,ind_out);                % C from gL and tau
            NeuPar(7,ind_out) = zeros(size(NeuPar(4,ind_out)));                      % a=0
            
            % detect parameters outside boundaries and draw them from a mulitvariate uniform distribution
            ind_out_para = [];
            for j=ind1
                ind_out_para = [ind_out_para find(NeuPar(ind3(j),iset)<N_min(ind3(j)+1,i) | NeuPar(ind3(j),iset) >N_max(ind3(j)+1,i))];  % check parameter boundaries
            end
            ind_out_V2 = find(NeuPar(9,iset)>=NeuPar(10,iset));    % Vr must be smaller than Vth
            ind_out_tau = find(NeuPar(1,iset)./NeuPar(2,iset)<N_min(12,i) | NeuPar(1,iset)./NeuPar(2,iset)>N_max(12,i)); % check tau boundaries
            ind_out_tcw = find(NeuPar(6,iset)<=NeuPar(1,iset)./NeuPar(2,iset)); % tcw must be larger than tau
            ind_out=iset(unique([ind_out_V2 ind_out_tau ind_out_tcw ind_out_para])); 
            iii=iii+1;
        end
        
        iii=0;
        while (~isempty(ind_out) && iii<1000)
            if iii==0
                disp('Number of trials for the full multivariate distribution has been exceeded. Use multivariate uniform distribution for the rest.');
            end
            % draw random numbers from a multivariate uniform distribution
            NeuPar_multi=multivar_distr(length(ind_out),N_min(ind3+1,i), N_max(ind3+1,i), ParCov{i},1)';
            
            % complete parameter matrix by "C", "a" and "Vup"
            for j=ind1
                NeuPar(ind3(j),ind_out) = NeuPar_multi(j,:);
            end
            NeuPar(1,ind_out) = NeuPar(1,ind_out).*NeuPar(2,ind_out);                % C from gL and tau
            NeuPar(7,ind_out) = zeros(size(NeuPar(4,ind_out)));                      % a=0
            
            % detect parameters outside boundaries and draw them from a mulitvariate uniform distribution
            ind_out_para = [];
            for j=ind1
                ind_out_para = [ind_out_para find(NeuPar(ind3(j),iset)<N_min(ind3(j)+1,i) | NeuPar(ind3(j),iset) >N_max(ind3(j)+1,i))];  % check parameter boundaries
            end
            ind_out_V2 = find(NeuPar(9,iset)>=NeuPar(10,iset));    % Vr must be smaller than Vth
            ind_out_tau = find(NeuPar(1,iset)./NeuPar(2,iset)<N_min(12,i) | NeuPar(1,iset)./NeuPar(2,iset)>N_max(12,i)); % check tau boundaries
            ind_out_tcw = find(NeuPar(6,iset)<=NeuPar(1,iset)./NeuPar(2,iset)); % tcw must be larger than tau
            ind_out=iset(unique([ind_out_V2 ind_out_tau ind_out_tcw ind_out_para])); 
            iii=iii+1;
        end
        
        if (~isempty(ind_out))
            error('Error: Number of trials for univariate uniform distribution has been exceeded.');
        end

%         for j = 2:length(SimPar.CellPar(:,1))      % use this block instead of the previous one before using update_inv_con_PSP
%             NeuPar(j-1,ind_out) = rand_par(length(ind_out), CellPar(j,i), N_sig(j,i), N_min(j,i), N_max(j,i), 0);
%         end;
%         NeuPar(3,[2:3 5:11 13:17])  = -60;                     % uniform EL above GABA reversal potential
%         NeuPar(5,[2:3 5:11 13:17])  = 10001; % no spiking in simulated neurons (except for the first one)
%         NeuPar(10,[2:3 5:11 13:17]) = 10000;

                   
        % Compute measures for electrophysiological classification of cells
        % and parameters for depolarization block 
        t_lat_act = zeros(size(iset));
        t_lat_LIF_act = zeros(size(iset));
        adapt_act = zeros(size(iset));
        for j=1:length(iset)

            % meaningful model parameter names
            [C,gL,EL,sf,~,~,~,~,Vr,Vth]=names(NeuPar(:,iset(j)));
            
            % compute parameters for depolarization block:
            % I_ref: Input current beyond which the neuron does not react
            %        any more
            % v_dep: Membrane potential that is used to compute the (input-
            %        independent) differential equation for the membrane
            %        potential during the depolarization block
            IRheo=gL*(Vth-EL)+gL*sf;
            OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
            I_ref = fminsearch(@Define_I_ref,IRheo+100,OPTIONS,NeuPar(:,iset(j)));
            NeuPar(11,iset(j)) = I_ref;     % I_ref
            NeuPar(12,iset(j)) = Vr;        % v_dep
            
            % compute latency of first spike
            I   = 500;
            t_lat_act(j) = integral(@(V) C./(I - gL*(V-EL) + gL*sf*exp((V-Vth)/sf) ), EL, Vth/2);
            t_lat_LIF_act(j) = C*log(I/(I+gL*(EL-Vth/2)))/gL;
%             t_lat_act(j) = 1;       % replace the two above lines by these ones before running update_inv_con_PSP
%             t_lat_LIF_act(j) = 1;
            
            
            % compute adaptation ratio
            I = 25:25:300;
            f1 = zeros(length(I),1);
            f  = zeros(length(I),1);
            for k=1:length(I)
                f1(k) = FRsimpAdEx(NeuPar(:,iset(j)),I(k),0,NeuPar(3,iset(j)),[]);
                f(k)  = FRsimpAdEx(NeuPar(:,iset(j)),I(k),[],[],[]); 
            end;
            if isnan(median(f1(f>0) ./ f(f>0)))
                adapt_act(j) = 0;
            else
                adapt_act(j) = median(f1(f>0) ./ f(f>0));
            end;
            
        end;
        t_lat(iset) = t_lat_act;
        t_lat_LIF(iset) = t_lat_LIF_act;
        adapt(iset) = adapt_act;
        
        % finally, set initial conditions
        V0(2,iset)=zeros(1,length(iset));
        for j=1:length(iset)
            V0(1,iset(j))=NeuPar(3,iset(j));   
        end
        
        NN = NN+NTypes(i);
    end;
    
    
    % Redistribute neuron types: I-L -> I-L-d             % comment out the redistribution before using update_inv_con_PSP
    ind_L23 = [reorder{2,ii} reorder{3,ii}];
    ind_L23_Ld = ind_L23((t_lat(ind_L23) - t_lat_LIF(ind_L23))>0);
    ind_L23_L = ind_L23(~ismember(ind_L23, ind_L23_Ld));
    NTypes(2) = length(ind_L23_L);
    NTypes(3) = length(ind_L23_Ld);
    reorder{2,ii} = ind_L23_L;
    reorder{3,ii} = ind_L23_Ld;
    
    ind_L5  = [reorder{9,ii} reorder{10,ii}];
    ind_L5_Ld = ind_L5((t_lat(ind_L5) - t_lat_LIF(ind_L5))>0);
    ind_L5_L = ind_L5(~ismember(ind_L5, ind_L5_Ld));
    NTypes(9) = length(ind_L5_L);
    NTypes(10) = length(ind_L5_Ld);
    reorder{9,ii} = ind_L5_L;
    reorder{10,ii} = ind_L5_Ld;
    
    
    % Redistribute neuron types: I-CL -> I-CL-AC
    ind_L23 = [reorder{4,ii} reorder{5,ii}];
    ind_L23_Ld = ind_L23(adapt(ind_L23)>1.5834);
    ind_L23_L = ind_L23(~ismember(ind_L23, ind_L23_Ld));
    NTypes(4) = length(ind_L23_L);
    NTypes(5) = length(ind_L23_Ld);
    reorder{4,ii} = ind_L23_L;
    reorder{5,ii} = ind_L23_Ld;
    
    ind_L5  = [reorder{11,ii} reorder{12,ii}];
    ind_L5_Ld = ind_L5(adapt(ind_L5)>1.5834);
    ind_L5_L = ind_L5(~ismember(ind_L5, ind_L5_Ld));
    NTypes(11) = length(ind_L5_L);
    NTypes(12) = length(ind_L5_Ld);
    reorder{11,ii} = ind_L5_L;
    reorder{12,ii} = ind_L5_Ld;

    

    % Set network connectivity
    if isempty(ConMtx)              % connect randomly if ConMtx is empty

        % intra-celltype connections
        for i=1:NTypesN
            if cluster_flag(i,i) == 1                                                                       % set random connectivity submatrix for the current input->output pair
                [X,idc] = SetCon_CommonNeighbour_Recur(Nsyn, length(reorder{i,ii}), length(reorder{i,ii}), pCon(i,i), 0.47, 1);  % with...
            else
                [X,idc] = SetCon(Nsyn, length(reorder{i,ii}), length(reorder{i,ii}), pCon(i,i));            % ... or without common neighbour rule
            end;
            
            for k = 1:length(ConPar{2}{i,i}{1})                                                             % loop over all synapse types
                l = ConPar{2}{i,i}{1}(k);                                                                   % determine which parameter set to use for each type
                ST = ConPar{2}{i,i}{2}(k);                                                                  % get synapse type ID (for SynPar)
                SPMtx(reorder{i,ii}, reorder{i,ii},k)=(X)+(k-1)*length(idc);                                % fit submatrix in
                
                idc_back = idc;
                [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,ConPar{2}{i,i}{l+2}, ...                        % set synapse parameters
                    ConPar{3}{i,i}, S_sig(i,i,:), S_max(i,i,:), S_min(i,i,:));
                if k==1
                    dtax_back = SynPar(6,idc_back);
                else
                    SynPar(6,idc_back) = dtax_back;
                    dtax_back = [];
                end;
            end;
        end;
        
        % inter-celltype connections
        for i=1:NTypesN                              % loop over output neurons
            for j=setdiff(1:NTypesN,i)               % loop over input neurons (except output neuron)
                if cluster_flag(i,j) == 2                                                                       % set random connectivity submatrix for the current input->output pair
                    [X,idc] = SetCon_CommonNeighbour_Cross(Nsyn, SPMtx(length(reorder{i,ii}),length(reorder{i,ii})), SPMtx(length(reorder{j,ii}),length(reorder{j,ii})), pCon(i,j), [], [], 1);  % with...
                    X=X';
                else
                    [X,idc] = SetCon(Nsyn, length(reorder{i,ii}), length(reorder{j,ii}), pCon(i,j));                                                                             % ... or without common neighbour rule
                end;
                
                for k = 1:length(ConPar{2}{i,j}{1})                                                             % loop over all synapse types
                    l = ConPar{2}{i,j}{1}(k);                                                                   % determine which parameter set to use for each type
                    ST = ConPar{2}{i,j}{2}(k);                                                                  % get synapse type ID (for SynPar)
                    SPMtx(reorder{i,ii}, reorder{j,ii},k)=(X)+(k-1)*length(idc);                                % fit submatrix in

                    idc_back = idc;
                        [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,ConPar{2}{i,j}{l+2}, ...                    % set synapse parameters
                            ConPar{3}{i,j}, S_sig(i,j,:), S_max(i,j,:), S_min(i,j,:));
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                end;
            end;
        end;

    else                            % use ConMtx if available
        
        for i=1:NTypesN                 % loop over input neurons
            for j=1:NTypesN             % loop over output neurons
                for k = 1:length(ConPar{2}{i,j}{1})                                     % loop over all (different) synapse types
                    X=ConMtx(reorder{i,ii}, reorder{j,ii},1);                           % Import connectivity submatrix for the current input->output pair
                    kk=find(X);                                                         % Extract non-zero entries...
                    idc=Nsyn+(1:length(kk));                                            % ... and create the idc vector from them
                    
                    X(kk)=idc;                                                          % Fill submatrix with idc entries...
                    X(X==0)=(k-1)*length(idc);                                          % ... and fill up the rest (for compatibitliy with former code only)
                    SPMtx(reorder{i,ii}, reorder{j,ii},k)=(X);                          % fit submatrix into full connectivity matrix
                    
                    l = ConPar{2}{i,j}{1}(k);                                           % determine which parameter set to use for each type
                    ST = ConPar{2}{i,j}{2}(k);                                          % get synapse type ID (for SynPar)

                    idc_back = idc;
                    [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,ConPar{2}{i,j}{l+2}, ...
                        ConPar{3}{i,j}, S_sig(i,j,:), S_max(i,j,:), S_min(i,j,:));      % set synapse parameters
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                end;
            end;
        end;
    end;

    
    % Add noise synapses (one for each synapse type)
    if max(max(NoisePar))>0
        for i=1:length(NoisePar(:,1))
            idc_back = idc;
            [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,length(SynPar(1,:))+1,i,NoisePar(i,:), {1.0,{[1 0 0], [0 0 0]}}, zeros(1,length(NoisePar(1,:))), zeros(1,length(NoisePar(1,:))),zeros(1,length(NoisePar(1,:))));
            if k==1
                dtax_back = SynPar(6,idc_back);
            else
                SynPar(6,idc_back) = dtax_back;
                dtax_back = [];
            end;
        end;
    end
end;


% ------------------------- Define inter-stripe connections ---------------------------------

for i=1:length(ConParStripes{1})
    i_act = ConParStripes{1}(i,1);
    j_act = ConParStripes{1}(i,2);
    
    for j=1:Nstripes
        for kk=1:length(ConParStripes{3}{i})
            
            % Set target stripe, dist_act and p_act
            target_act = j + ConParStripes{3}{i}(kk);
            while target_act<=0
                target_act = target_act + Nstripes;
            end;
            while target_act>Nstripes
                target_act = target_act - Nstripes;
            end;
            dist_act = abs(ConParStripes{3}{i}(kk));
            p_act = pCon(i_act,j_act)*exp(-dist_act/ConParStripes{2}{i}(1));

            
            % insert matrix and parameter definitions (as for intra-stripe connections)
            if isempty(ConMtx)              % connect randomly if ConMtx is empty
                if cluster_flag(i_act,j_act) == 2                                                                       % set random connectivity submatrix for the current input->output pair
                    [X,idc] = SetCon_CommonNeighbour_Cross(Nsyn, SPMtx(length(reorder{i_act,target_act}),length(reorder{i_act,target_act})), SPMtx(length(reorder{j_act,j}),length(reorder{j_act,j})), p_act, 1);   % with...
                    X=X';
                else
                    [X,idc] = SetCon(Nsyn, length(reorder{i_act,target_act}), length(reorder{j_act,j}), p_act);                                                                                                     % ... or without common neighbour rule
                end;
                
                for k = 1:length(ConPar{2}{i_act,j_act}{1})                                                             % loop over all synapse types
                    l = ConPar{2}{i_act,j_act}{1}(k);                                                                   % determine which parameter set to use for each type
                    ST = ConPar{2}{i_act,j_act}{2}(k);                                                                  % get synapse type ID (for SynPar)
                    SPMtx(reorder{i_act,target_act}, reorder{j_act,j},k)=(X)+(k-1)*length(idc);                         % fit submatrix in
                    
                    % Define synapse parameters from horizontal distance
                    wgt_act = ConPar{2}{i_act,j_act}{l+2}(1)*exp(-dist_act/ConParStripes{2}{i}(1)); 
                    dtax_act = ConPar{2}{i_act,j_act}{l+2}(2) + ConParStripes{2}{i}(2)*dist_act;
                    S_sig_act(1) = S_sig(i_act,j_act,1)*exp(-dist_act/ConParStripes{2}{i}(1));
                    S_sig_act(2) = S_sig(i_act,j_act,2) + ConParStripes{2}{i}(2)*dist_act;
                    
                    idc_back = idc;
                    [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,[wgt_act dtax_act], ...                                 % set synapse parameters
                        ConPar{3}{i_act,j_act}, S_sig_act, S_max(i_act,j_act,:), S_min(i_act,j_act,:));
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                end;
            else                            % use ConMtx if available

                for k = 1:length(ConPar{2}{i_act,j_act}{1})                                                 % loop over all (different) synapse types
                    X=ConMtx(reorder{i_act,target_act}, reorder{j_act,j},1);                                % Import connectivity submatrix for the current input->output pair
                    kkk=find(X);                                                                            % Extract non-zero entries...
                    idc=Nsyn+(1:length(kkk));                                                               % ... and create the idc vector from them
                    
                    X(kkk)=idc;                                                                             % Fill submatrix with idc entries...
                    X(X==0)=(k-1)*length(idc);                                                              % ... and fill up the rest (for compatibitliy with former code only)
                    SPMtx(reorder{i_act,target_act}, reorder{i_act,j},k)=(X);                               % fit submatrix into full connectivity matrix
                    
                    l = ConPar{2}{i_act,j_act}{1}(k);                                                       % determine which parameter set to use for each type
                    ST = ConPar{2}{i_act,j_act}{2}(k);                                                      % get synapse type ID (for SynPar)

                    % Define synapse parameters from horizontal distance
                    wgt_act = ConPar{2}{i_act,j_act}{l+2}(1)*exp(-dist_act/ConParStripes{2}{i}(1)); 
                    dtax_act = ConPar{2}{i_act,j_act}{l+2}(2) + ConParStripes{2}{i}(2)*dist_act;
                    S_sig_act(1) = S_sig(i_act,j_act,1)*exp(-dist_act/ConParStripes{2}{i}(1));
                    S_sig_act(2) = S_sig(i_act,j_act,2) + ConParStripes{2}{i}(2)*dist_act;
                    
                    idc_back = idc;
                    [SynPar,Nsyn,~]=SetSyn(SynPar,Nsyn,idc,ST,[wgt_act dtax_act], ...
                        ConPar{3}{i_act,j_act}, S_sig_act, S_max(i_act,j_act,:), S_min(i_act,j_act,:));      % set synapse parameters
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                end;
            end;
        end;
    end;
end;


% ------------------------- Define connections to input neurons ---------------------------------

for ii=1:Nstripes
    reorder{NTypesN+1,ii} = NN+1:NN+InpNTypes(1);
    reorder{NTypesN+2,ii} = NN+InpNTypes(1)+1:NN+M;
    
    if isempty(ConMtx)              % connect randomly if ConMtx is empty
        
        % inter-celltype connections
        for i=1:NTypesN                              % loop over output neurons
            for j=NTypesN+1:NTypesN+InpNTypesN       % loop over input neurons (literally)
                if cluster_flag(i,j) == 2                                                                       % set random connectivity submatrix for the current input->output pair
                    [X,idc] = SetCon_CommonNeighbour_Cross(Nsyn, SPMtx(length(reorder{i,ii}),length(reorder{i,ii})), SPMtx(length(reorder{j,ii}),length(reorder{j,ii})), pCon(i,j), [], [], 1);  % with...
                    X=X';
                else
                    [X,idc] = SetCon(Nsyn, length(reorder{i,ii}), length(reorder{j,ii}), pCon(i,j));                                                                             % ... or without common neighbour rule
                end;
                
                for k = 1:length(ConPar{2}{i,j}{1})                                                             % loop over all synapse types
                    l = ConPar{2}{i,j}{1}(k);                                                                   % determine which parameter set to use for each type
                    ST = ConPar{2}{i,j}{2}(k);                                                                  % get synapse type ID (for SynPar)
                    SPMtx(reorder{i,ii}, reorder{j,ii},k)=(X)+(k-1)*length(idc);                                % fit submatrix in
                    
                    idc_back = idc;
                    [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,ConPar{2}{i,j}{l+2}, ...                        % set synapse parameters
                        ConPar{3}{i,j}, S_sig(i,j,:), S_max(i,j,:), S_min(i,j,:));
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                    
                end;
            end;
        end;
        
    else                            % use ConMtx if available
        
        for i=1:NTypesN                              % loop over input neurons
            for j=NTypesN+1:NTypesN+InpNTypesN       % loop over input neurons (literally)
                for k = 1:length(ConPar{2}{i,j}{1})                                     % loop over all (different) synapse types
                    X=ConMtx(reorder{i,ii}, reorder{j,ii},1);                           % Import connectivity submatrix for the current input->output pair
                    kk=find(X);                                                         % Extract non-zero entries...
                    idc=Nsyn+(1:length(kk));                                            % ... and create the idc vector from them
                    
                    X(kk)=idc;                                                          % Fill submatrix with idc entries...
                    X(X==0)=(k-1)*length(idc);                                          % ... and fill up the rest (for compatibitliy with former code only)
                    SPMtx(reorder{i,ii}, reorder{j,ii},k)=(X);                          % fit submatrix into full connectivity matrix
                    
                    l = ConPar{2}{i,j}{1}(k);                                           % determine which parameter set to use for each type
                    ST = ConPar{2}{i,j}{2}(k);                                          % get synapse type ID (for SynPar)
                    
                    idc_back = idc;
                    [SynPar,Nsyn,~]=SetSyn(SynPar,Nsyn,idc,ST,ConPar{2}{i,j}{l+2}, ...
                        ConPar{3}{i,j}, S_sig(i,j,:), S_max(i,j,:), S_min(i,j,:));      % set synapse parameters
                    if k==1
                        dtax_back = SynPar(6,idc_back);
                    else
                        SynPar(6,idc_back) = dtax_back;
                        dtax_back = [];
                    end;
                end;
            end;
        end;
    end;
end



% ----------------------------- Define cell assemblies --------------------------------------
if ~isempty(CA)
    SPMtxCA = SPMtx(CA(:,1),CA(:,1),:);

    for i=1:length(CA(:,1))
        CA_ind = find(SPMtxCA(:,i,1)>0);
        if(~isempty(CA_ind))
            for k=1:2       % exc. synapses only
                SynPar(5,SPMtxCA(CA_ind,i,k)) = CA(i,2)*SynPar(5,SPMtxCA(CA_ind,i,k));
            end;
        end;
    end;
end;


% ------------------------- Save full parameter set  ---------------------------------
if ~isempty(SimPar.fnOut)
    EvtMtx=SimPar.EvtMtx;
    EvtTimes=SimPar.EvtTimes;
    save(SimPar.fnOut,'SimPar', 'CtrPar','NeuPar','NPList', ...
        'STypPar','SynPar','SPMtx','EvtMtx','EvtTimes','ViewList','V0', '-v7.3');
end;


% -------------- Run the actual simulation (described in IDNet.c) -------------------------

[Ninp,~]=IDNet(CtrPar,NeuPar,NPList,STypPar,SynPar,SPMtx,SimPar.EvtMtx,SimPar.EvtTimes,ViewList,InpSTtrains,NoiseDistr,V0,UniqueNum,NeuronGroupsSaveArray);


% -------------- Convert data into binary MATLAB format -------------------------
% check if there are neurons where the input oscillates around I_ref (not necessary any more)
% idx_osc=find(N_osc>0);
% if ~isempty(idx_osc)
%     disp('The following neurons might have a problem with their input (oscillates around I_ref):');
%     disp(num2str(idx_osc));
% end

% Read in temporary files and delete them
STMtx=[];
if CtrPar(5)>0
    for i=1:N*Nstripes+M
        fnST=['ISIu' num2str(i-1) '_' num2str(UniqueNum) '.dat'];
        if exist(fnST,'file')
            ST=load(fnST);
            STMtx{i}=ST';
            delete(fnST);
        else
            STMtx{i}=[];
        end;
    end;
end;

T=[]; V=[];dV=[];
if ~isempty(ViewList)
    X=load(['IDN_' num2str(UniqueNum) '.dat']);
    nv=length(ViewList);
    T=X(1:nv:end,1);
    for i=1:nv
        V{i}=X(i:nv:end,3:end);
    end;
    delete(['IDN_' num2str(UniqueNum) '.dat'])
end;

if exist(['IDN2_' num2str(UniqueNum) '.dat'],'file')
    dX=load(['IDN2_' num2str(UniqueNum) '.dat']);
    for Counter1 = 1:3*length(NeuronGroupsSaveArray)
        dV{Counter1} = dX(:,Counter1);
    end
    delete(['IDN2_' num2str(UniqueNum) '.dat'])
end

if exist(['IDN3_' num2str(UniqueNum) '.dat'],'file')
    XX=importdata(['IDN3_' num2str(UniqueNum) '.dat']);
    delete(['IDN3_' num2str(UniqueNum) '.dat'])
end


% ------------------------- Save results  ---------------------------------
if ~isempty(SimPar.fnOut)
    EvtMtx=SimPar.EvtMtx;
    EvtTimes=SimPar.EvtTimes;
    save(SimPar.fnOut,'Ninp','T','V','dV','XX','STMtx', '-append', '-v7.3');
end;



% -------------------------------------------------------------------------
% -------------- Functions to generate connectivity -----------------------
% -------------------------------------------------------------------------

function [X,idc]=SetCon(Nsyn, Nin,Nout,pCon)
% Generates a random connectivity matrix X[ConPar(1), ConPar(2)] with
% connection prob. ConPar(3). Nsyn (number of synapses) is continuated
% from input matrices. idc gives the running indices of the generated
% connections (reshaped in one dimension)
% 
% INPUT:
%   Nsyn:   Number of synapses that are already created
%   Nin:    Number of input neurons
%   Nout:   Number of output neurons
%   pcon:   Connection probability
% 
% OUTPUT:
%   X:      Connection matrix (Nin x Nout matrix) containing idc on 
%           non-zero entries
%   idc:    Synapse indices (starting with Nsyn)


k1 = randperm(Nin*Nout);
k = k1(1:round(pCon*Nin*Nout));
X=zeros(Nin,Nout);
idc=Nsyn+(1:length(k));
X(k)=idc;


function [SynPar,Nsyn,idc]=SetSyn(SynPar,Nsyn,idc,ST,ConPar,ConParSTSP,S_sig,S_max,S_min)
% Sets (i.e. continuates) synaptic parameter matrix SynPar(6,idc), where
% idc is the running index of connections (see SetCon). 
% 
% INPUT:
%   SynPar:         Synapse parameters (excluding those computed here)
%   Nsyn:           Number of synapses that are already set
%   idc:            Indices of synapses that are already set
%   ST:             Synapse type
%   ConPar:         Mean connection parameters (synaptic weight, delay and 
%                   failure probability)
%   ConParSTSP:     Short-term synaptic plasticity parameters
%   S_sig:          Standard deviation of connection parameters
%   S_max:          Maximum of connection parameters
%   S_min:          Minimum of connection parameters
% 
% OUTPUT:
%   SynPar:         Synapse parameters (concatinating old and new columns)
%   Nsyn:           Updated number of synapses
%   idc:            Updated synapse indices

wgt=ConPar(1)';
dtax=ConPar(2);
SynPar(1,idc) = ST;                                                        % Synapse type (AMPA, GABA and NMDA here)


mean_wgt = log((wgt^2)/sqrt(S_sig(1)^2+wgt^2));
std_wgt = sqrt(log(S_sig(1)^2/(wgt^2)+1));
if wgt==0 || S_sig(1)==0
    SynPar(5,idc) = rand_par(length(idc), wgt, S_sig(1),   S_min(1)*wgt,     S_max(1)*wgt,    0);            % wgt
else
    SynPar(5,idc) = rand_par(length(idc), mean_wgt,    std_wgt,    S_min(1)*wgt,     S_max(1)*wgt,    2);    % wgt
end

SynPar(6,idc) = rand_par(length(idc), dtax,        S_sig(2),   S_min(2)*dtax,    S_max(2)*dtax,   0);        % dtax
SynPar(7,idc) = ConPar(3);                                                                                   % p_fail  

if ~isempty(idc)
    % Deal with rounding errors
    idc_frac = round(ConParSTSP{1}*length(idc));
    if sum(idc_frac) > length(idc)
        [~,ind] = min(ConParSTSP{1}*length(idc) - floor(ConParSTSP{1}*length(idc)));
        idc_frac(ind) = idc_frac(ind)-1;
    else
        if sum(idc_frac) < length(idc)
            [~,ind] = max(idc_frac - floor(idc_frac));
            idc_frac(ind) = idc_frac(ind)+1;
        end;
    end;
    
    idc_act = Nsyn+(1:idc_frac(1));
    Nsyn_act = Nsyn+idc_frac(1);
    for i=1:length(ConParSTSP{1})
        use    = ConParSTSP{i+1}{1}(1);
        tc_rec = ConParSTSP{i+1}{1}(2);
        tc_fac = ConParSTSP{i+1}{1}(3);
        std_use    = ConParSTSP{i+1}{2}(1);
        std_tc_rec = ConParSTSP{i+1}{2}(2);
        std_tc_fac = ConParSTSP{i+1}{2}(3);
        
        % set maxima and standard deviations to zero before using update_inv_con_PSP
        SynPar(2,idc_act) = rand_par(length(idc_act), use,    std_use,    0,    1, 0);    % use
        SynPar(3,idc_act) = rand_par(length(idc_act), tc_rec, std_tc_rec, 0, 1500, 0);    % tc_rec
        SynPar(4,idc_act) = rand_par(length(idc_act), tc_fac, std_tc_fac, 0, 1500, 0);    % tc_fac
        
        if i<length(ConParSTSP{1}) && length(ConParSTSP{1})>1
            idc_act=Nsyn_act+(1:idc_frac(i+1));
            Nsyn_act=Nsyn_act+idc_frac(i+1);
        end;
    end;
end;

Nsyn=Nsyn+length(idc);
idc=Nsyn+(1:length(idc));


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
