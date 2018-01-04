function q=Define_I_ref(I0,Par)
% Objective function that needs to be minimized to compute the current I
% above which the cell does not react any more (depolarization block)
%
% INPUT:
%    I0:    Current in pA
%    Par:   Model parameters
%
% OUTPUT
%    q:     The value of the objective function

    re=FRsimpAdEx(Par,I0,0,[],[]);
    q = (re-200).^2;

end


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
