% Create hazard curve
% Density function of the selected conditional(?) variable

function hazard_pdf = hazard(Control, Structure)

scale           = Control.fragi.condi_var.scale_factor;
var_index       = strncmp(Control.fragi.condi_var.name, Structure.probdata.name, 100);
condi_var_pos   = var_index == 1;

marg            = Structure.probdata.marg;
marg            = marg(condi_var_pos,:);
dist_type       = marg(1);

% Gumbel!!
%     u_n             = marg(5)*scale; % WARNING
%     a_n             = marg(6)/scale; % WARNING
%     hazard_pdf      = @(x)a_n * exp( -a_n*(x-u_n) - exp(-a_n*(x-u_n))); %Gumbel max! pdf

% more general solution using FERUM
% it can happen that initialization of marg is missing.. this partially corrects that
if isnan(marg(8))
    marg(8) = 1;
end
% WARNING! scaling is added!
hazard_pdf      = @(x) ferum_pdf(dist_type, x/scale , marg(5:8))/scale;

end