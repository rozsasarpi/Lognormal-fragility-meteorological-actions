% Fit requested distribution to a sample using least square regression in untransformed space (LS)
%
%SYNOPSIS
% [cdf, x] = FIT_LS(fragi_sample, intensity, distr)
%
%INPUT
%
%OUTPUT
% cdf           cumulative distribution function, /numeric vector, nx1/
% x             associated values, /numeric vector, nx1/
%
%NOTE(S)
% * The range of cdf is based on the intensity vector of reliability based fragility curve
%   that already should be siufficiently wide and fine, but here it more refined and exteneded for sure
% * The empirical distribution function is not unique, herein
%   empiF_i = (i - 2/5)/(k + 1/5) [Cunnane, 1978]
% SEE ALSO
% fit_ML, fit_MM, fit_LST

function [cdf, x] = fit_LS(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        
        obs     = fragi_sample;
        sobs    = sort(obs);
        k       = length(obs);

        empiF   = (((1:k) - 2/5)/(k + 1/5))';
        
        % method of moments estimate for initial value
        mu_logn         = log(median(fragi_sample));
        sigma_logn      = sqrt(log((std(fragi_sample)).^2./mean(fragi_sample).^2 + 1));
        x0              = [mu_logn, sigma_logn];

        cdf_fun         = @(x, xdata) logncdf(xdata, x(1), x(2));
        
%         if any(isnan(x0)) || any(isinf(x0))
%             error('Nye!')
%         end
        options = optimoptions('lsqcurvefit', 'Display', 'off');
        % least square fit
        [parm_hat, ~, ~, exitflag]  = lsqcurvefit(cdf_fun, x0, sobs, empiF, [], [], options);
        
        if exitflag < 1
            warning('Problem with the numerical optimization!')
        end
        
        cdf             = logncdf(x, parm_hat(1), parm_hat(2));
        
%         % DIAGNOSTIC PLOT
%         plot(sobs, empiF, 'o')
%         hold on
%         plot(x, cdf)
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end