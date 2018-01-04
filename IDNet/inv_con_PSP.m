function gmax = inv_con_PSP(PSP,i,j)
% Implements linear transformation from a desired maximal PSP value to a peak
% conductance, according to the data in inv_con_par, depending on input and
% output neuron type

% scaling factors for all target neuron types
par_E = [1.0569    0.5875    0.6587    0.7567    0.6728    0.9899    0.6294    1.6596    0.5941    0.6661    0.7647    0.6799    1.5818    0.6360];
par_I = [2.3859    1.6277    1.6277    1.6671    1.6671    2.3142    1.4363    3.5816    1.6277    1.6277    1.6671    1.6671    3.4016    1.4363];
% par_E = ones(1,14);  % replace the former two lines with these ones before using update_inv_con_PSP
% par_I = ones(1,14);


if any(ismember([1 8 15],j))       % use par_E if input is excitatory
    gmax = PSP*par_E(i);
else                               % use par_I if input is inhibitory
    gmax = PSP*par_I(i); 
end;
