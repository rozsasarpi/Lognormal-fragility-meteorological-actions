% Check the effect of lognormal assumption for fragility curves
% fit models and produce publication quality plots
% quite general implementation
%
% ra - reliability analysis based
%

%NOTES:
% Variable naming convention:
%  Structure - captial first letter, lower case the rest
%  MATRIX    - matrix (dim>=2) all caps (not always..)
%
% Structure - levels hierarchy: from larger sets (outer) to smaller ones (inner/deeper)
% results are available as multidimensional 
%
clear variables
close all
clc

%--------------------------------------------------------------------------
% MAIN CONTROL PARAMETERS
%--------------------------------------------------------------------------
n_boot = 1e3;   % number of bootstrap simulations to estimate sampling variability
n_simu = 42;    % number of samples/simulations to estimate fragility curve parameters (lognormal)
p_ci   = 0.9;   % confidence interval level

% only lognormal is available yet, but the code is prepared with extensibiltiy in mind
distr  = 'ln';  % distribution type for fragility function

% RA:   Original result of reliable analysis
% F0:   Simple fitted result whitout sampling
% MM:   method of moments (Porter Method A)
% MLM:  method of L-moments
% ML:   maximum likelihood method
% TML:  truncated sample + ML
% LS:   least square regression
% LST:  least square regression in transformed space (Porter Method B)
fit_methods = {'F0', 'MM', 'MLM', 'LS', 'LST', 'ML'};
% fit_methods = {'F0', 'MM', 'MLM', 'LST', 'ML'};

set(0,'defaulttextinterpreter','latex')
cmp = get(groot,'defaultAxesColorOrder');

%--------------------------------------------------------------------------
% LOAD DATA & INITIALIZE VARIABLES
%--------------------------------------------------------------------------
% turbine hall - system
% load('D:\Working folder\Matlab working folder\Nubiki fragility\Számítások dokumentálása\Acél\keretállás_rendszer\Turbinacsarnok\analysis_IO_teljes_rendszer.mat')

% reactor hall - system
% load('J:\MELOO\PAKS\MATHLAB_szamitasok\REAKTOR_CSARNOK\__RENDSZER_SZINT\FRAGI_KET_KOMP_KERETALLAS\analysis_IO+FRAGI_KET KOMPLETT KERETALLASRA+Keret_6rendszer_lsf_V1_olok_MIN1_szel_E_intervall_ge_szil_rendszer_DOKSI_UJ_SZEL.mat')
% load('E:\DOKTORI\CIKKEK_ELOADASOK\Fragility_cikk_2016\automatizalando_file/analysis_IO_teljes_rendszer.mat')
% load([pwd, '\analysis_IO_teljes_rendszer.mat'])
load('d:\Working folder\Matlab working folder\SS_lognormal assumption\turbinacsarnok\data\cikk_16_32a_uj.mat')

% fRAgi outputs
% WARNING!
if exist('Control2', 'var')
    Control     = Control2;
    Structure   = Structure2;
    Results     = Results2;
end
% WARNING!

type        = 'full';
% beta_ra     = Results.fragi.(type).system.muB_bounds(:,1);
beta_ra     = Results.fragi.(type).system.beta_exact(:,1);
% beta_ra     = Results.fragi.(type).g(3).beta;
cdf_ra      = normcdf(-beta_ra);
scale       = Control.fragi.condi_var.scale_factor;
intensity   = Results.fragi.intensity*scale;

% Populate structure variable
% common for all fragiltiy curves/ fitting methods
Fit.Common.intensity    = intensity;
Fit.Common.scale        = scale;
Fit.Common.p_ci         = p_ci;

Fit.RA.fragi_cdf        = cdf_ra;

%--------------------------------------------------------------------------
% DEFINE THE DISTRIBUTION TYPE
%--------------------------------------------------------------------------
switch lower(distr)
    case {'ln', 'ln2', 'lnorm', 'lognorm', 'lognormal'}
        type_cdf = @(x,par1,par2)(logncdf(x,par1,par2));
    otherwise
        error(['Unknown distribution ID: ', distr])
end


%--------------------------------------------------------------------------
% GET MOMENTS AND OTHER CHARACTERISTICS
%--------------------------------------------------------------------------

% overwrites _cdf_ra_ and _intensity_ with a 'cleaned' vesion: no duplicates in cdf and pdf
Fit = get_statistics(Fit);

% populate the workspace with one-level deep structure fields
cellfun(@(field) assignin('caller', field, Fit.Common.(field)), fieldnames(Fit.Common))

%--------------------------------------------------------------------------
% FITTING METHODS & INITIALIZATION
%--------------------------------------------------------------------------

all_methods = ['RA',fit_methods];
n_method = length(fit_methods);

% extended and refine intensity range for integration
eintensity      = extend_intensity(intensity);
n_i             = length(eintensity);
Fit.Common.eintensity  = eintensity;

% uniformly distributed random numbers for bootstrap
rng(1983)
ru              = rand(n_simu, n_boot);


%--------------------------------------------------------------------------
% HAZARD CURVE & INTEGRATION FOR RA
%--------------------------------------------------------------------------

% from script: % WARNING scale, see script, or els... but the hazard_pdf
% must be an anonym function
hazard_pdf_fun  = hazard(Control,Structure);
hazard_pdf      = hazard_pdf_fun(intensity);

% calculate failure probability
fragi_cdf   = Fit.RA.fragi_cdf;
dpf         = hazard_pdf.*fragi_cdf;
pf_ra       = trapz(intensity, dpf);
Fit.RA.dpf  = dpf;
Fit.RA.pf   = pf_ra;
Fit.RA.beta = -norminv(pf_ra);

% --------------------------------------------------------------------------
% 'FIT' THE CHOOSEN DISTRIBUTION TYPE BY DIFFERENT METHODS
% --------------------------------------------------------------------------

%Initialization of bootsrap
% dim1 (rows)       intensity
% dim2 (columns)    bootstrap sample
% dim3 (depth)      method
FRAGI_CDF       = nan(n_i, n_boot, n_method);
X               = nan(n_i, n_boot, n_method);
DPF             = nan(n_i, n_boot, n_method);
PF              = nan(1, n_boot, n_method);
PF_RATIO        = nan(1, n_boot, n_method);
hazard_pdf      = hazard_pdf_fun(eintensity);

% loop over fitting methods
% it unnecessarily bootstraps F0 as well
for kk = 1:n_method
    fit_method = fit_methods{kk};
    fhandle = ['fit_',fit_method];
    
    % bootstrap to estimate sampling variability
    % loop over bootstrap sample
    parfor ii = 1:n_boot
        
        % fit distribution
        % WARNING!
%         fragi_sample = interp1(cdf_ra, intensity, ru(:,ii), 'linear','extrap');
        fragi_sample = interp1(cdf_ra, intensity, ru(:,ii));
%         if any(isnan(fragi_sample))
%             error('Nye!')
%         end
        if kk == 1
            [fragi_cdf, x] = feval(fhandle, mean_ra, median_ra, std_ra, intensity, distr);
        else
            [fragi_cdf, x] = feval(fhandle, fragi_sample, intensity, distr);
        end
        
        % calculate failure probability
        dpf         = hazard_pdf.*fragi_cdf;
        pf          = trapz(eintensity, dpf);
        pf_ratio    = pf/pf_ra;
        
        FRAGI_CDF(:,ii,kk)   = fragi_cdf;
        X(:,ii,kk)           = x;
        DPF(:,ii,kk)         = dpf;
        PF(:,ii,kk)          = pf;
        PF_RATIO(:,ii,kk)    = pf_ratio;
    end
    
    Fit.(fit_method).fragi_cdf  = FRAGI_CDF(:,:,kk);
    
    Fit.(fit_method).dpf        = DPF(:,:,kk);
    Fit.(fit_method).pf         = PF(:,:,kk);
    Fit.(fit_method).beta       = -norminv(PF(:,:,kk));
    Fit.(fit_method).pf_ratio   = PF_RATIO(:,:,kk);
    Fit.(fit_method).ci         = quantile(FRAGI_CDF(:,:,kk), [(1-p_ci)/2, (1+p_ci)/2], 2);
    Fit.(fit_method).median     = median(FRAGI_CDF(:,:,kk), 2);
    Fit.(fit_method).mean       = mean(FRAGI_CDF(:,:,kk), 2);
    
end

%--------------------------------------------------------------------------
% PLOT - SMALL MULTIPLES - OPTIONS
%--------------------------------------------------------------------------

% initializing subplot
% Options:
Options.fontsize    = 9;
Options.imp_range   = 1; % imp: importance range 0:No 1:Yes
Options.US          = 1; % US: untransformed space; 0:No 1:Yes
Options.LNS         = 1; % LNS: lognormal space; 0:No 1:Yes
Options.F0          = 1; % F0: Simple 'type' fit 0:No 1:Yes
Options.subsize     = 3; % number of figure next to each other (2 or 3 is possible)
Options.xlim        = [0,15];


%--------------------------------------------------------------------------
% PLOT - SMALL MULTIPLES - UNTRANSFORMED
%--------------------------------------------------------------------------

if Options.US == 1
    % initializing subplot
    [Pos, plot_methods] = subplot_fit(fit_methods, Options);
    plot_fragi(Fit,Options, Pos, plot_methods);
end

%=======================================
% SAVE
%=======================================
% if saveplot == 1
%     % get style sheet info
%     snam   = 'ss_smallm_3.4';    % name of style file (NO extension)
%     s      = hgexport('readstyle',snam);
%     
%     %apply style sheet info
%     fnam   = ['d:\Working folder\Matlab working folder\SS_lognormal assumption\turbinacsarnok\figures\_', tag,'.eps']; % file name
%     s.Format = 'eps';
%     hgexport(gcf, fnam, s);
% else
%     % do nothing
% end
% ------------------------------------------------
% PLOT THE LOGNORMAL FIGURE
% ------------------------------------------------

if Options.LNS == 1
    % initializing subplot
    [Pos,plot_methods] = subplot_fit_LN(fit_methods, Options);
    plot_fragi_LN(Fit,Options ,Pos, plot_methods);
end


