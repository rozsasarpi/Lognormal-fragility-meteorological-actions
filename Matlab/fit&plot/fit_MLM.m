% Fit requested distribution to a sample using method of L-moments (MLM)
%
%SYNOPSIS
% [cdf, x] = FIT_MLM(fragi_sample, intensity, distr)
%
%INPUT
%
%OUTPUT
% cdf           cumulative distribution function, /numeric vector, nx1/
% x             associated values, /numeric vector, nx1/
%
%NOTE
% the range of cdf is based on the intensity vector of reliability based fragility curve
% that already should be siufficiently wide and fine, but here it more refined and exteneded for sure
%
% SEE ALSO
% fit_ML, fit_MM, fit_LS, fit_LST

function [cdf, x] = fit_MLM(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        % L-moments of the sample
        L               = lmom(fragi_sample, 2);
        sigma_logn      = 2*erfinv(L(2)/L(1));          % [Singh, 1998, Eq.(6.61)]
        mu_logn         = log(L(1)) - sigma_logn^2/2;   % [Singh, 1998, Eq.(6.62)]
        
        
        cdf             = logncdf(x, mu_logn, sigma_logn);  
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end