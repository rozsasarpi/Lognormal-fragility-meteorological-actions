% Check the effect of lognormal assumption for Paks steel hall(s)
%
% ra - reliability analysis based
%
% first crude implementation; to be wrapped into small functions

clear all
close all
clc

N_boot = 1e3; % number of bootstrap simulations to estimate sampling variability
n_simu = 44;  % number of samples/simulations to estimate fragility curve parameters (lognormal)
p_ci   = 0.9; % confidence interval level

set(0,'defaulttextinterpreter','latex')
cmp = get(groot,'defaultAxesColorOrder');

%--------------------------------------------------------------------------
% LOAD DATA
%--------------------------------------------------------------------------
% turbine hall - system
load('D:\Working folder\Matlab working folder\Nubiki fragility\Számítások dokumentálása\Acél\keretállás_rendszer\Turbinacsarnok\analysis_IO_teljes_rendszer.mat')

% reactor hall - system
% load('J:\MELOO\PAKS\MATHLAB_szamitasok\REAKTOR_CSARNOK\__RENDSZER_SZINT\FRAGI_KET_KOMP_KERETALLAS\analysis_IO+FRAGI_KET KOMPLETT KERETALLASRA+Keret_6rendszer_lsf_V1_olok_MIN1_szel_E_intervall_ge_szil_rendszer_DOKSI_UJ_SZEL.mat')
% load('E:\DOKTORI\CIKKEK_ELOADASOK\Fragility_cikk_2016\automatizalando_file/analysis_IO_teljes_rendszer.mat')

% WARNING!
Control = Control2;
Structure = Structure2;
Results = Results2;
% WARNING!

beta_ra     = Results.fragi.separated.system.muB_bounds(:,1);
P_ra        = normcdf(-beta_ra);
scale       = Control.fragi.condi_var.scale_factor;
intensity   = Results.fragi.intensity*scale;

%--------------------------------------------------------------------------
% GET MOMENTS AND OTHER CHARACTERISTICS
%--------------------------------------------------------------------------
fragi_cdf_ra    = @(x) interp1(intensity, P_ra, x);
fragi_pdf_ra    = @(x) cfd(fragi_cdf_ra, x);
p_ra            = fragi_pdf_ra(intensity);
p_ra(isnan(p_ra)) = 0;

% masking the same value
[P_ra,ia]       = unique(P_ra);
intensity       = intensity(ia);
p_ra            = p_ra(ia);

mean_ra         = trapz(intensity, intensity.*p_ra);
median_ra       = interp1(P_ra, intensity, 0.5);
std_ra          = trapz(intensity, (mean_ra - intensity).^2.*p_ra);

%--------------------------------------------------------------------------
% 'FIT' LOGNORMAL CURVES
%--------------------------------------------------------------------------
%.........................................................................
% A1
mu_logn         = log(median_ra);
sigma_logn      = sqrt(log((std_ra).^2./mean_ra.^2 + 1));
% fragi_cdf_a1    = @(x) lognormcdf(x, mean_ra, std_ra/mean_ra);
fragi_cdf_a1    = @(x) logncdf(x, mu_logn, sigma_logn);

parmhat_B1 = nan(2, N_boot);
parmhat_B2 = nan(2, N_boot);
% bootstrap to estimate sampling variability
parfor ii = 1:N_boot
    %.........................................................................
    % B1 - method of moments
    rng(332+ii) % for reproducibility
    ru = rand(n_simu, 1);
    fragi_sample = interp1(P_ra, intensity, ru);
    
    mu_logn         = log(median(fragi_sample));
    sigma_logn      = sqrt(log((std(fragi_sample)).^2./mean(fragi_sample).^2 + 1));
    % fragi_cdf_b1    = @(x) logncdf(x, mu_logn, sigma_logn);
    
    parmhat_B1(:,ii) = [mu_logn, sigma_logn];
    %..........................................................................
    % B2 - max likelihood
    parmhat         = fit_lognorm2_mle(fragi_sample, 'par');
    % fragi_cdf_b2    = @(x) logncdf(x, parmhat(1), parmhat(2));
    
    parmhat_B2(:,ii) = parmhat;
    
end
%--------------------------------------------------------------------------
% INTEGRATION
%--------------------------------------------------------------------------
% using the firs bootsrtap element
fragi_cdf_b1    = @(x) logncdf(x, parmhat_B1(1,1), parmhat_B1(2,1));
fragi_cdf_b2    = @(x) logncdf(x, parmhat_B2(1,1), parmhat_B2(2,1));

var_index       = strncmp(Control.fragi.condi_var.name, Structure.probdata.name, 100);
condi_var_pos   = find(var_index == 1);

marg            = Structure.probdata.marg;
dist_type       = marg(condi_var_pos,1);

u_n             = marg(condi_var_pos,5)*scale; % WARNING
a_n             = marg(condi_var_pos,6)/scale; % WARNING

% Gumbel!!
hazard_pdf      = @(x)a_n * exp( -a_n*(x-u_n) - exp(-a_n*(x-u_n))); %Gumbel max! pdf
d_Pf_ra         = @(x)hazard_pdf(x).*fragi_cdf_ra(x);
d_Pf_a1         = @(x)hazard_pdf(x).*fragi_cdf_a1(x);
d_Pf_b1         = @(x)hazard_pdf(x).*fragi_cdf_b1(x);
d_Pf_b2         = @(x)hazard_pdf(x).*fragi_cdf_b2(x);

Pf_ra           = integral(d_Pf_ra, intensity(1), intensity(end));     % WARNING - depends on the chosen range of the intensity!
Pf_a1           = integral(d_Pf_a1, intensity(1), intensity(end));     % WARNING - depends on the chosen range of the intensity!
Pf_b1           = integral(d_Pf_b1, intensity(1), intensity(end));     % WARNING - depends on the chosen range of the intensity!
Pf_b2           = integral(d_Pf_b2, intensity(1), intensity(end));     % WARNING - depends on the chosen range of the intensity!


Beta_ra         = norminv(1-Pf_ra);
Beta_a1         = norminv(1-Pf_a1);
Beta_b1         = norminv(1-Pf_b1);
Beta_b2         = norminv(1-Pf_b2);


[Pf_a1, median(Pf_b1), median(Pf_b2)]./Pf_ra

%--------------------------------------------------------------------------
% PLOT - SMALL MULTIPLES - UNTRANSFORMED
%--------------------------------------------------------------------------
figure('Position', [500, 500, 800, 400])
%..........................................................................
% A1
%..........................................................................
subplot(2,3,1)
h_ra = plot(intensity, P_ra, 'black');
hold on
h_a1 = plot(intensity, fragi_cdf_a1(intensity), 'black--');

text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{A1}} = ', num2str(roundsd(Beta_a1,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('A1')
ylabel('$P_\mathrm{f}(q)$')
xlabel('Ground snow $s\,[kN/m^2]$')
xlim([0,15])

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')

%..........................................................................
% B1
%..........................................................................
% boostrap confidence intervals
P_b1        = nan(N_boot, length(intensity));
Pf_ratio_b1 = nan(N_boot,1);
parfor ii = 1:N_boot
    fragi           = @(x) logncdf(x, parmhat_B1(1,ii), parmhat_B1(2,ii));
    P_b1(ii,:)      = fragi(intensity);
    d_Pf_b1         = @(x)hazard_pdf(x).*fragi(x);
    Pf_b1           = integral(d_Pf_b1, intensity(1), intensity(end));
    Pf_ratio_b1(ii) = Pf_b1/Pf_ra;
end
ci_b1 = quantile(P_b1, [(1-p_ci)/2, (1+p_ci)/2])';
b_med1= median(P_b1)';

subplot(2,3,2)
h_ra = plot(intensity, P_ra, 'black');
hold on
h_a1 = plot(intensity, fragi_cdf_b1(intensity), 'black--');
[~, hb] = boundedline(intensity, b_med1, [b_med1 - ci_b1(:,1),  ci_b1(:,2) - b_med1]);
uistack(hb, 'bottom')

text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{B1}} = ', num2str(roundsd(Beta_b1,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('B1')
% ylabel('$P_\mathrm{f}(q)$')
set(gca,'YTickLabel','')
xlabel('Ground snow $s\,[kN/m^2]$')
xlim([0,15])
ylim([0,1])

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')


% Failure probability ratio plot - TODO: fill under the curve with color of confidence bounds
subplot(2,3,5)
if median(Pf_ratio_b1) < 1
    ss = 1./Pf_ratio_b1;
    tag = '$P_\mathrm{f,RA}/P_\mathrm{f,B1}$';
else
    ss = Pf_ratio_b1;
    tag = '$P_\mathrm{f,B1}/P_\mathrm{f,RA}$';
end
[~, density, x] = kde(ss);
plot(x, density, 'black')
xlabel(tag)
ylabel('density')
xlim([0,100])
hold on
mea = mean(ss);
med = median(ss);
plot([mea, mea], ylim, 'black')
plot([med, med], ylim, 'black')
ht = text(mea-0.05*diff(xlim), mean(ylim)-0.2*diff(ylim), ['$mean = ', num2str(roundsd(mea,3)), '$']);
ht.Rotation = 90;
ht = text(med-0.05*diff(xlim), mean(ylim)-0.2*diff(ylim), ['$median = ', num2str(roundsd(med,3)), '$']);
ht.Rotation = 90;

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')


%..........................................................................
% B2
%..........................................................................
% boostrap confidence intervals
P_b2        = nan(N_boot, length(intensity));
Pf_ratio_b2 = nan(N_boot,1);
parfor ii = 1:N_boot
    fragi           = @(x) logncdf(x, parmhat_B2(1,ii), parmhat_B2(2,ii));
    P_b2(ii,:)      = fragi(intensity);
    d_Pf_b2         = @(x)hazard_pdf(x).*fragi(x);
    Pf_b2           = integral(d_Pf_b2, intensity(1), intensity(end));
    Pf_ratio_b2(ii) = Pf_b2/Pf_ra;
end
ci_b2 = quantile(P_b2, [(1-p_ci)/2, (1+p_ci)/2])';
b_med2 = median(P_b2)';

subplot(2,3,3)
h_ra = plot(intensity, P_ra, 'black');
hold on
h_a1 = plot(intensity, fragi_cdf_b2(intensity), 'black--');
[~, hb] = boundedline(intensity, b_med2, [b_med2 - ci_b2(:,1),  ci_b2(:,2) - b_med2]);
uistack(hb, 'bottom')

text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{B2}} = ', num2str(roundsd(Beta_b2,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('B2')
% ylabel('$P_\mathrm{f}(q)$')
set(gca,'YTickLabel','')
xlabel('Ground snow $s\,[kN/m^2]$')
xlim([0,15])
ylim([0,1])

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')


% Failure probability ratio plot - TODO: fill under the curve with color of confidence bounds
subplot(2,3,6)
if median(Pf_ratio_b2) < 1
    ss = 1./Pf_ratio_b2;
    tag = '$P_\mathrm{f,RA}/P_\mathrm{f,B2}$';
else
    ss = Pf_ratio_b2;
    tag = '$P_\mathrm{f,B2}/P_\mathrm{f,RA}$';
end
[~, density, x] = kde(ss);
plot(x, density, 'black')
xlabel(tag)
ylabel('density')
xlim([0,100])
hold on
mea = mean(ss);
med = median(ss);
plot([mea, mea], ylim, 'black')
plot([med, med], ylim, 'black')
ht = text(mea-0.05*diff(xlim), mean(ylim)-0.2*diff(ylim), ['$mean = ', num2str(roundsd(mea,3)), '$']);
ht.Rotation = 90;
ht = text(med-0.05*diff(xlim), mean(ylim)-0.2*diff(ylim), ['$median = ', num2str(roundsd(med,3)), '$']);
ht.Rotation = 90;

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')

%--------------------------------------------------------------------------
% PLOT - SMALL MULTIPLES - TRANSFORMED
%--------------------------------------------------------------------------
xx      = log(intensity);

figure('Position', [500, 500, 800, 400])
%..........................................................................
% A1
%..........................................................................
subplot(2,3,1)
h_ra = plot(xx, erfinv(2*P_ra-1), 'black');
hold on
h_a1 = plot(xx, erfinv(2*fragi_cdf_a1(intensity)-1), 'black--');

ylim([-6,6])
ylim_a1 = get(gca,'Ylim');

text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{A1}} = ', num2str(roundsd(Beta_a1,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('A1')
ylabel('$erf^{-1}\left(2P_\mathrm{f}(q)-1 \right)$')
xlabel('Ground snow $log(s\,[kN/m^2])$')

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')

%..........................................................................
% B1
%..........................................................................
subplot(2,3,2)
h_ra    = plot(xx, erfinv(2*P_ra-1), 'black');
hold on
h_a1    = plot(xx, erfinv(2*fragi_cdf_b1(intensity)-1), 'black--');
tci_b1  = erfinv(2*ci_b1-1);
tb_med1 = erfinv(2*b_med1-1);

idx     = logical(prod(~isinf(tci_b1),2));
[~, hb] = boundedline(xx(idx), tb_med1(idx), [tb_med1(idx) - tci_b1(idx,1),  tci_b1(idx,2) - tb_med1(idx)]);
uistack(hb, 'bottom')

ylim(ylim_a1)
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{B1}} = ', num2str(roundsd(Beta_b1,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('B1')
% ylabel('$erf^{-1}\left(2P_\mathrm{f}(q)-1 \right)$')
set(gca,'YTickLabel','')

xlabel('Ground snow $log(s\,[kN/m^2])$')    

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')

%..........................................................................
% B2
%..........................................................................
subplot(2,3,3)
h_ra = plot(xx, erfinv(2*P_ra-1), 'black');
hold on
h_a1 = plot(xx, erfinv(2*fragi_cdf_b2(intensity)-1), 'black--');
tci_b2  = erfinv(2*ci_b2-1);
tb_med2 = erfinv(2*b_med2-1);

idx     = logical(prod(~isinf(tci_b2),2));
[~, hb] = boundedline(xx(idx), tb_med2(idx), [tb_med2(idx) - tci_b2(idx,1),  tci_b2(idx,2) - tb_med2(idx)]);
uistack(hb, 'bottom')

ylim(ylim_a1)
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
    ['$\beta_{\mathrm{RA}} = ', num2str(roundsd(Beta_ra,3)),'$'])
text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
    ['$\beta_{\mathrm{B2}} = ', num2str(roundsd(Beta_b2,3)),'$'])

% legend([h_ra, h_a1], {'RA', 'A1'})

title('B2')
% ylabel('$erf^{-1}\left(2P_\mathrm{f}(q)-1 \right)$')
set(gca,'YTickLabel','')
ylim(ylim_a1)
xlabel('Ground snow $log(s\,[kN/m^2])$')

set(gca,'TickDir','out')
% set(gca,'TickLabelInterpreter', 'LaTeX')