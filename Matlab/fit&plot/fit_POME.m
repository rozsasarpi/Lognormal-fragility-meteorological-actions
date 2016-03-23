% Fit requested distribution to a sample using maximum entropy principle with fixed distribution type (POME)
%
%SYNOPSIS
% [cdf, x] = FIT_POME(fragi_sample, intensity, distr)
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


function [cdf, x] = fit_POME(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        % [Singh, 1998, 6.1.5. chapter]
        % same as the maximum likelihood estimate, I am not sure, second look is needed..
% %         n           = numel(fragi_sample);
% %         mu_hat      = sum(log(fragi_sample))/n;
% %         sigma_hat   = sqrt(sum(log(fragi_sample).^2)/n - mu_hat^2);
% % 
% %         cdf         = logncdf(x, mu_hat, sigma_hat);  
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end