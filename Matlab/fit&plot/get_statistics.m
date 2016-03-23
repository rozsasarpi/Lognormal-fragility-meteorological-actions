% Get statistics of a random variable given with cumulative distribution function 
%
%SYNOPSIS
% Fit = GET_STATISTICS(x, cdf)
%
%INPUT
% x             /numeric vector, nx1/
% cdf           cumulative distribution function /numeric vector, nx1/
% 
%OUTPUT
% Fit.
%  Common.
%   mean_ra 
%   median_ra
%   std_ra
%   pdf_ra
%   cdf_ra
%   inensity
%
% SEE ALSO

function Fit = get_statistics(Fit)

cdf               = Fit.RA.fragi_cdf;
x                 = Fit.Common.intensity;

fragi_cdf_ra      = @(x0) interp1(x, cdf, x0);
fragi_pdf_ra      = @(x0) cfd(fragi_cdf_ra, x0);
pdf_ra            = fragi_pdf_ra(x);
pdf_ra(isnan(pdf_ra)) = 0;

% masking the same values
[cdf_ra, ia]     = unique(cdf);
x_ra             = x(ia);
pdf_ra           = pdf_ra(ia);

mean_ra          = trapz(x_ra, x_ra.*pdf_ra);
median_ra        = interp1(cdf_ra, x_ra, 0.5);
std_ra           = sqrt(trapz(x_ra, (mean_ra - x_ra).^2.*pdf_ra));


% Define strucutre of variables
Fit.Common.mean_ra     = mean_ra;
Fit.Common.median_ra   = median_ra;
Fit.Common.std_ra      = std_ra;
Fit.Common.pdf_ra      = pdf_ra;
Fit.Common.cdf_ra      = cdf_ra;
Fit.Common.intensity   = x_ra;

Fit.RA.fragi_cdf       = cdf_ra;

end