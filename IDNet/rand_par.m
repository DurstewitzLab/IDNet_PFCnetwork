function par = rand_par(N, par_mean, par_std, par_min, par_max, distr_flag)
% Returns N random parameters. If a drawn parameter falls outside the 
% interval [par_min par_max], it is replaced by a random variable drawn 
% from a uniform distribution spanning that interval. % If par_std = 0 and 
% par_max = 0, all parameters are set to par_mean.

% 
% INPUT:
%   N:          Number of random parameters
%   par_mean:   Parameter mean
%   par_std:    Parameter standard deviation
%   par_min:    Parameter minimum
%   par_max:    Parameter maximum
%   distr_flag: Specifies the distribution used for the random draw:
%                   0: Gaussian distribution
%                   1: uniform distribution over interval [par_min par_max]
%                   2: log-normal distribution
% 
% OUTPUT:
%   par:        Vector of random parameters


if par_std==0
    if par_max==0
        par = par_mean*ones(1,N);                   % constant parameters with std and par_max == 0
    else
        par = par_min+(par_max-par_min)*rand(1,N);  % uniform distribution if std==0, but par_max !=0
    end;
else
    switch distr_flag                               % choose distribution according to distr_flag
        case 0
            par = par_mean + par_std*randn(1,N);
            exc_ind = find(par<par_min | par>par_max);
            par(exc_ind) =  par_min+(par_max-par_min)*rand(1,length(exc_ind));
        case 1
            par = par_min+(par_max-par_min)*rand(1,N);
        case 2
            par = exp(randn(N,1) .* par_std + par_mean);
            exc_ind = find(par<par_min | par>par_max);
            par(exc_ind) =  par_min+(par_max-par_min)*rand(1,length(exc_ind));            
    end;
end;

% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
