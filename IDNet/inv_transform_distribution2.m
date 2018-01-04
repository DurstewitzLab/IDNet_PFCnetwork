function X_inv=inv_transform_distribution2(X_trans,k,mean_X,std_X,min_X)
% inverts the transformation made in transform_distribution
% using the moments of the original distribution
%
% INPUT:
%   X_trans:    Sample from the transformed distribution
%   k:          Box-Cox exponent
%   mean_X:     Mean of the transformed distribution
%   std_X:      Standard deviation of the transformed distribution
%   min_X:      Minimum of the transformed distribution
% 
% OUTPUT:
%   X_inv:      Sample from the inverted distribution


std_decr = 0.8;

if k>0
    if min_X<0
        X_inv = (mean_X+std_decr*std_X*X_trans).^(1/k)+1.1*min_X;
    else
        X_inv = (mean_X+std_decr*std_X*X_trans).^(1/k);
    end
else
    if min_X<0
        X_inv = exp(mean_X+std_decr*std_X*X_trans)+1.1*min_X;
    else
        X_inv = exp(mean_X+std_decr*std_X*X_trans);
    end
end


% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
