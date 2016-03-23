% Fit requested distribution to a sample using least square regression in lognormal space (LST)
%
%SYNOPSIS
% [cdf, x] = FIT_LST(fragi_sample, intensity, distr)
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
% fit_ML, fit_MM, fit_LS

function [cdf, x] = fit_LST(fragi_sample, intensity, distr)

x = extend_intensity(intensity);

switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        
        obs             = fragi_sample;
        sobs            = sort(obs);
        k               = length(obs);

        empiF           = (((1:k) - 2/5)/(k + 1/5))';
        
        % transform to lognormal space
        tsobs           = log(sobs);
        tempiF          = erfinv(2*empiF-1);
        
        % least square fit
        p               = polyfit(tempiF, tsobs, 1);
        parm_hat        = [p(2), p(1)/sqrt(2)];
        
        cdf             = logncdf(x, parm_hat(1), parm_hat(2));
        
%         % DIAGNOSTIC PLOT
%         plot(tsobs, tempiF, 'o')
%         hold on
%         plot(log(x), erfinv(2*cdf-1))
    otherwise
        error(['Unknown distribution ID: ', distr])
end

end