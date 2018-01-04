function [par_E, par_I23, par_I5, res_E, res_I23, res_I5] = update_inv_con_PSP(gmax)
% update the connections strengths that lead to a defined PSP
% amplitude, using a range of gmax values (optional input)
%
% INPUT:
%   gmax: Peak conductances to be tested (vector)
% 
% OUTPUT:
%   p_X:   Parameters of the linear fit to the gmax/PSP function 
%          (X=E: excitatory input, X=I23: inhibitory input, L2/3, X=L5:
%          inhibitory input, L5)
%   res_X: gmax and PSP values, X as above


% Set test parameter set with only one neuron
SimParTest = ConfigIDNet(14,1,1,1000,0,[1.0 1.0],[1.0 1.0]);                   % only one neuron per type

SimParTest.EvtMtx = zeros(17,1);                    % input only to three neurons (pyramidal cell and interneuon in both layers)
SimParTest.EvtMtx(1,1) = 58.5;
SimParTest.EvtMtx(4,1) = 31.5;
SimParTest.EvtMtx(12,1) = 31.5;
SimParTest.EvtTimes = [SimParTest.CtrPar(3); SimParTest.CtrPar(2)];

SimParTest.N_sig = zeros(size(SimParTest.N_sig));   % no variation in parameters
SimParTest.N_min = SimParTest.CellPar;
SimParTest.N_max = zeros(size(SimParTest.N_max));
SimParTest.S_sig = zeros(size(SimParTest.S_sig));
SimParTest.S_max = zeros(size(SimParTest.S_max));


% Default values of gmax
if ~exist('gmax','var')
    gmax = 0:0.25:2;
end;

% Generate PSP data for all values of gmax (excitatory input)
SimParTest.ViewList = [2:3 5:11 13:17];                          % watch all neurons except the spiking ones
for i=1:length(gmax)
    for j=1:14
        SimParTest.ConPar{2}{j,1}{3}(1) = gmax(i);               % input strength
    end;
    M = zeros(17,19);                                            % Define connectivity from spiking neuron
    M(SimParTest.ViewList,1) = 1;                                % to all others
    [~,T,V,~,~]=IDNetSim(SimPar,M);
    res_E(i,1) = gmax(i);
    for j=1:length(SimParTest.ViewList)
        res_E(i,j+1) = max(V{j}(T>300)) - mean(V{j}(T>300 & T<400));
    end;
    disp([i res_E(i,:)])
end;

% Linear fit
for i=1:length(res_E(1,:))-1
    par_E(:,i) = polyfit(res_E(:,i+1), res_E(:,1), 1)';
end;


% Generate PSP data for all values of gmax (inhibitory input, L2/3)
SimParTest.ViewList = [2:3 5:9];                                 % watch all neurons in L2/3 except the spiking ones
for i=1:length(gmax)
    for j=1:14
        SimParTest.ConPar{2}{j,3}{3}(1) = gmax(i);               % input strength
    end;
    M = zeros(17,19);  % [0 1 0]
    M(SimParTest.ViewList,4) = 1;
    [~,T,V,~,~]=IDNetSim(SimPar,M);
    res_I23(i,1) = gmax(i);
    for j=1:length(SimParTest.ViewList)
        res_I23(i,j+1) = mean(V{j}(T>300 & T<400)) - min(V{j}(T>300));
    end;
    dips([i res_I23(i,:)])
end;

% Linear fit
for i=1:length(res_I23(1,:))-1
    par_I23(:,i) = polyfit(res_I23(:,i+1), res_I23(:,1), 1)';
end;


% Generate PSP data for all values of gmax (inhibitory input, L5)
SimParTest.ViewList = [10:11 13:17];                              % watch all neurons in L5 except the spiking ones
for i=1:length(gmax)
    for j=1:14
        SimParTest.ConPar{2}{j,10}{3}(1) = gmax(i);               % input strength
    end;
    M = zeros(17,19);  % [0 1 0]
    M(SimParTest.ViewList,12) = 1;
    [~,T,V,~,~]=IDNetSim(SimPar,M);
    res_I5(i,1) = gmax(i);
    for j=1:length(SimParTest.ViewList)
        res_I5(i,j+1) = mean(V{j}(T>300 & T<400)) - min(V{j}(T>300));
    end;
    disp([i res_I5(i,:)])
end;

% Linear fit
for i=1:length(res_I5(1,:))-1
    par_I5(:,i) = polyfit(res_I5(:,i+1), res_I5(:,1), 1)';
end;


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
