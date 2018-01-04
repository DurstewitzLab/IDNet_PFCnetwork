%% Compute new scaling factors for inv_con_PSP using update_inv_con_PSP 
%
% Run this whenever parameters are changed that affect the post-synaptic
% potential (PSP) in any of the neuron classes.
% Use par_E and par_I computed here in inv_con_PSP.m
%
%
% IMPORTANT:
% Before running this script, several changes need to be made:
%
%
% update_inv_con_PSP (rarely):
% - Adjust SimParTest.EvtMtx(1,1) such that the first neuron elicits a
%   single spike at around 500 ms (only if neuron parameters where changed).
%
% In test_update_inv_con_PSP (rarely):
% - Change STP values if they were changed in ConfigIDNet
%
% In ConfigIDNet (always):
% - Change NTypes line in ConfigIDNet to 
%    "NTypes = ones(1,14); NTypes(1) = 2;"
%
% In IDNetSim (always):
% - Prevent multivariate random number drawings by the line  
%    "if ismember(i,[])" 
% - Comment out the block starting with 
%    "while ~isempty(ind_out)" 
%   and use the previous block instead
% - Set "t_lat_act(j) = 1" and t_lat_LIF_act(j) = 1
%    instead of the real calculation
% - Comment out redistribution of neuron types
% - Set maxima and stds of STP parameters to zero in SetSyn
% - Set mean STP values to E2 or to I2, respectively, wherever there is more 
%   than one option
%
% In inv_con_PSP (always):
% - Set all par_E and par_I to one


% Compute scaling factors from simulation
gmax = [0.1 0.5 1.0 1.5 2.5 3.0];           % gmax to be tested
[par_E_act, par_I23_act, par_I5_act, res_E, res_I23, res_I5] = update_inv_con_PSP(gmax);

par_E = [(1:14)' par_E_act(1,:)'];
par_I = [(1:14)' [par_I23_act(1,:)'; par_I5_act(1,:)']];


% Gauge scaling factors with relative fractions of STP types (use only data
% within layers, since STP data across layers is only avaiable for E->E
% connections, which are almost exclusively depressing)
% D=Depressing, F=Facilitating, DF=Mixed
load('STP_types.mat')

mean_E1 = 0.28;   % Change these values if they are changed in ConfigIDNet
mean_E2 = 0.25;
mean_E3 = 0.29;
mean_I1 = 0.16;
mean_I2 = 0.25;
mean_I3 = 0.32;

par_E_scaled(1)  = (D_L23E_L23E + mean_E1/mean_E2*F_L23E_L23E + mean_E3/mean_E2*DF_L23E_L23E)*par_E(1,2);  % E->L23E

par_E_scaled(2)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(2,2);  % E->L23I-L
par_E_scaled(3)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(3,2);  % E->L23I-L-d
par_E_scaled(4)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(4,2);  % E->L23I-CL
par_E_scaled(5)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(5,2);  % E->L23I-CL-AC
par_E_scaled(6)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(6,2);  % E->L23I-CS
par_E_scaled(7)  = (D_L23E_L23I + mean_E1/mean_E2*F_L23E_L23I + mean_E3/mean_E2*DF_L23E_L23I)*par_E(7,2);  % E->L23I-F


par_E_scaled(8)  = (D_L5E_L5E + mean_E1/mean_E2*F_L5E_L5E + mean_E3/mean_E2*DF_L5E_L5E)*par_E(8,2);  % E->L5E

par_E_scaled(9)  = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(9,2);  % E->L5I-L
par_E_scaled(10) = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(10,2); % E->L5I-L-d
par_E_scaled(11) = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(11,2); % E->L5I-CL
par_E_scaled(12) = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(12,2); % E->L5I-CL-AC
par_E_scaled(13) = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(13,2); % E->L5I-CS
par_E_scaled(14) = (D_L5E_L5I + mean_E1/mean_E2*F_L5E_L5I + mean_E3/mean_E2*DF_L5E_L5I)*par_E(14,2); % E->L5I-F


par_I_scaled(1)  = (D_L23I_L23E + mean_E1/mean_E2*F_L23I_L23E + mean_E3/mean_E2*DF_L23I_L23E)*par_I(1,2);  % I->L23E

par_I_scaled(2)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(2,2);  % I->L23I-L
par_I_scaled(3)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(3,2);  % I->L23I-L-d
par_I_scaled(4)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(4,2);  % I->L23I-CL
par_I_scaled(5)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(5,2);  % I->L23I-CL-AC
par_I_scaled(6)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(6,2);  % I->L23I-CS
par_I_scaled(7)  = (D_L23I_L23I + mean_E1/mean_E2*F_L23I_L23I + mean_E3/mean_E2*DF_L23I_L23I)*par_I(7,2);  % I->L23I-F

par_I_scaled(8)  = (D_L5I_L5E + mean_E1/mean_E2*F_L5I_L5E + mean_E3/mean_E2*DF_L5I_L5E)*par_I(8,2);  % I->L5E

par_I_scaled(9)  = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(9,2);  % I->L5I-L
par_I_scaled(10) = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(10,2); % I->L5I-L-d
par_I_scaled(11) = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(11,2); % I->L5I-CL
par_I_scaled(12) = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(12,2); % I->L5I-CL-AC
par_I_scaled(13) = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(13,2); % I->L5I-CS
par_I_scaled(14) = (D_L5I_L5I + mean_E1/mean_E2*F_L5I_L5I + mean_E3/mean_E2*DF_L5I_L5I)*par_I(14,2); % I->L5I-F


par_E = [par_E par_E_scaled']; disp(par_E);
par_I = [par_I par_I_scaled']; disp(par_I);

save inv_con_PSP_factors par_E par_I


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
