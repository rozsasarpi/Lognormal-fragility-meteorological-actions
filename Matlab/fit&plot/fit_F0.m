% Fit requested distribution by matching moments exactly
%
%SYNOPSIS
% [cdf, x] = FIT_F0(mean_ra, median_ra, std_ra, intensity, distr)
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

function [cdf, x] = fit_F0(mean_ra, median_ra, std_ra, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        % match median
        mu_logn     = log(median_ra);
        % match mean
%         mu_logn     = log(mean_ra.^2./sqrt(std_ra.^2 + mean_ra.^2));
        
        sigma_logn  = sqrt(log((std_ra).^2./mean_ra.^2 + 1));
        
        cdf         = logncdf(x, mu_logn, sigma_logn);  
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end