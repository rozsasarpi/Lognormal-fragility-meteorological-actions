% Fit requested distribution to a sample using maximum likelihood method (ML)
%
%SYNOPSIS
% [cdf, x] = FIT_ML(fragi_sample, intensity, distr)
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
% fit_lognorm2_mle


function [cdf, x] = fit_ML(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        parm_hat    = fit_lognorm2_mle(fragi_sample, 'par');
        
        cdf         = logncdf(x, parm_hat(1), parm_hat(2));  
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end