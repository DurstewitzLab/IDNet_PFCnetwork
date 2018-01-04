%% Worksheet runing the simulations and creating a raster plot  of the 
%% resulting spike train

% Main output: 
% - STMtx: Spike times in ms, (cell array, one cell per neuron)
% - V: Membrane potential in mV, (cell array, one cell per neuron)
% - T: Simulation time in ms (vector)


% -----------------  1) Prepare simulations  ------------------

% Set paths
addpath('IDNet');

% Compile MEX file (only needed at first run, and when IDNet.c is changed)
cd('IDNet')
mex IDNet.c
cd('..')

% Set overall simulation parameters
N=1000;         % Number of neurons
M = 1;          % Number of input neurons
Str=1;          % Number of columns/stripes
SimTim=1000;    % Simulation time in ms
T_skip=500;     % Initial part of the spike train to skip for analysis in ms

sEE=1; sIE=1; sEI=1; sII=1;                         % Synaptic weight scales
pEE = 1; pEI = 1; pIE = 1; pII = 1; pE=1; pI=1;     % Connectivity scales
I = zeros(1,14);                                    % Background input 
I(1) = 250;                                         % I_ex
I(8) = 250;
I(2:7) = 200;                                       % I_inh
I(9:14) = 200;

% Compute complete simulation parameter set and construct file name
SimPar = ConfigIDNet(N,M,Str,SimTim,I,[sEE sIE sEI sII],[pEE pEI pIE pII pE pI]);
filename_1=['PFC_' num2str(I(1)) '_' num2str(I(2)) '_' num2str(N) 'N_S' num2str(Str)];
filename_2=['_s_' num2str(sEE*1) '_' num2str(sIE*1) '_' num2str(sEI*1) '_' num2str(sII*1)];
filename_3=['_p_' num2str(pEE) '_' num2str(pEI) '_' num2str(pIE) '_' num2str(pII) '_' num2str(pE) '_' num2str(pI) '_' num2str(SimTim) 'ms'];
filename=[filename_1 filename_2 filename_3];


% ---------  2)  Perform simulation or load existing data  ------------
if ~exist(filename,'file') && ~exist([filename '_all.mat'], 'file')
    tic;
    disp(filename)
    SimPar.fnOut=[filename '_all'];
    SimPar.CtrPar(2) = SimTim;
    SimPar.ViewList = 1;%1:sum(SimPar.NTypes);          % Set of neurons to record the membrane potential from
    SimPar.NeuronGroupsSaveArray=[];%SimPar.ViewList;   % Set to view list to monitor additional variables
    [STMtx,T,V,dV,Ninp]=IDNetSim(SimPar,[]);
    t=toc;
    disp(['Simulation time:' num2str(t) 's'])
else
    load([filename '_all'])
end


% --------  3) Display spike trains  ---------
STMtx=STMtx(1,1:length(STMtx)-2);
STMtx=cellfun(@(x) x(x>T_skip),STMtx,'UniformOutput',false);

figure;
for ii=1:length(STMtx)
    hold on
    plot(STMtx{ii}, ii*ones(length(STMtx{ii}),1), 'b.')
end;


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
