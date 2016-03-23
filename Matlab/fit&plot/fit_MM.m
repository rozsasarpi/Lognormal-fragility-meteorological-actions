% Fit requested distribution to a sample using method of moments (MM)
%
%SYNOPSIS
% [cdf, x] = FIT_MM(fragi_sample, intensity, distr)
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
% fit_ML

function [cdf, x] = fit_MM(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        mean_fragi      = mean(fragi_sample);
        median_fragi    = median(fragi_sample);
        std_fragi       = std(fragi_sample);
        
        % match median
%         mu_logn         = log(median_fragi);
        % match mean
        mu_logn         = log(mean_fragi.^2./sqrt(std_fragi.^2 + mean_fragi.^2));
        
        sigma_logn      = sqrt(log(std_fragi.^2./mean_fragi.^2 + 1));
        
        cdf             = logncdf(x, mu_logn, sigma_logn);  
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end