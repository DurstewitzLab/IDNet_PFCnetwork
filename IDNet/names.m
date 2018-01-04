function [Cm,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]=names(para)
% Assign meaningful names to all parameters
%
% INPUT:
%   para:   Vector of parameters
% 
% OUTPUT:
%   Cm:     Membrane capacity
%   gL:     Leak conductance
%   EL:     Lead reversal potential
%   sf:     Slope factor of the exponential part
%   Vup:    Maximal membrane potential (hard threshold)
%   tcw:    Adaptation time constant
%   a:      Continuous adaptation parameter 
%   b:      Spike-triggered adaptation parameter
%   Vr:     Reversal potential
%   Vth:    Threshold potential (soft threshold - exponential term sets in)


    Cm=para(1);
    gL=para(2);
    EL=para(3);
    sf=para(4);
    Vup=para(5);
    tcw=para(6);
    a=para(7);
    b=para(8);
    Vr=para(9);
    Vth=para(10);
    
    
% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
