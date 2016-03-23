% Plot figures in lognormal space
%

function plot_fragi_LN(Fit, Options, Pos, plot_methods)

% ------------------------------------------------
% PLOT THE LOGNORMAL FIGURE
% ------------------------------------------------

n_methods = length(plot_methods);
n_col = Options.subsize;
n_row = ceil(n_methods/n_col);

figure('Position', [400, 200, n_col*200, n_row*250])

for ii = 1:n_methods
    
    plot_method = plot_methods{ii};
    intensity   = Fit.Common.intensity;
    eintensity  = Fit.Common.eintensity;
%     scale      = Fit.Common.scale;
    p_ci       = Fit.Common.p_ci;
    cdf_ra     = Fit.Common.cdf_ra;
    fit_cdf    = Fit.(plot_method).fragi_cdf(:,1);

    logi        = log(intensity);
    logei       = log(eintensity);
    yy_ra       = erfinv(2*cdf_ra-1);
    yy_fit_met  = erfinv(2*fit_cdf-1);
    if ~strcmp(plot_method,'F0')
%         conf_med   = erfinv(2*Fit.conf.(plot_method).med-1);
%         conf_ci    = erfinv(2*Fit.conf.(plot_method).ci-1);
        % to trick erfinv
        conf_med   = norminv(Fit.(plot_method).median, 0, 1/sqrt(2));
        conf_ci    = norminv(Fit.(plot_method).ci, 0, 1/sqrt(2));
    end
    
    % fragility curve
    axes('position',Pos{ii},'XGrid','off','XMinorGrid','off','FontSize',Options.fontsize,'Box','on','Layer','top');
    
        %..........................................................................
        % IMPORTANCE CURVE
        %..........................................................................
        if Options.imp_range == 1
                
            importance   = Fit.RA.dpf/Fit.RA.pf;

            importance_cdf = cumtrapz(intensity, importance);
            [~, ia, ~] = unique(importance_cdf);
            ci_lower = interp1(importance_cdf(ia), intensity(ia), (1-p_ci)/2);
            ci_upper = interp1(importance_cdf(ia), intensity(ia), (1+p_ci)/2);
            
            ci_lower = log(ci_lower);
            ci_upper = log(ci_upper);
            
            
            hold on
            patch([ci_lower, ci_lower, ci_upper, ci_upper], [-5, 5, 5, -5], zeros(4,1),...
                                            'EdgeColor', 'none',...
                                            'FaceColor', [1, 0.8, 0.8])
            set(gca, 'Layer','top')

        end
        
        % Add bounds
        if ~strcmp(plot_method,'F0')
            idx     = logical(prod(~isinf(conf_ci),2));            
            [hl, hb] = boundedline(logei(idx), conf_med(idx), [conf_med(idx) - conf_ci(idx,1),  conf_ci(idx,2) - conf_med(idx)]);
            uistack([hl,hb], 'top')
            
        end
   
        plot(logi,yy_ra, 'black');
        hold on
        plot(logei, yy_fit_met, 'black--');

 
        ylim([-4,5])

        title(plot_method)
        xlabel('Ground snow $log(s)\,[kN/m^2]$')
        
        xlm = Options.xlim(2);
        digit = 10; %10 means 0.1; 100 means 0.01
        
        xlim(log([round(digit*xlm/50)/digit,Options.xlim(2)]))

        set(gca,'TickDir','out')
        set(gca,'TickLabelInterpreter', 'LaTeX')
        % grid on
        
        % Ticks
        x_tick = [round(digit*xlm/50)/digit, 1, 2, 3, 5, Options.xlim(2)];
%       x_tick      = get(gca,'XTick');
        x_tick_transf   = log(x_tick);
        set(gca,'XTick',x_tick_transf);
        set(gca,'XTickLabel',x_tick );
        xlim(log([round(digit*xlm/50)/digit,Options.xlim(2)]))

%         y_tick          = [10^-15, 10^-10, 10^-5, 10^-3, 10^-2, 10^-1, 0.5, 0.9, 0.999, 0.999999];
        y_tick          = [10^-15, 10^-10, 10^-5, 10^-3, 10^-1, 0.5, 0.9, 0.999];
        if Pos{ii}(1) <= 0.20
            ylabel('$erf^{-1}\left(2P_\mathrm{f}(s)-1 \right)$')
            y_tick_transf   = erfinv(2*y_tick-1);
            set(gca,'YTick',y_tick_transf)
%             set(gca,'YTickLabel', y_tick)
            set(gca,'YTickLabel', sprintf('%0.2e\n',y_tick))
        else
            y_tick_transf   = erfinv(2*y_tick-1);
            set(gca,'YTick',y_tick_transf)
            set(gca,'YTickLabel','')    
        end

        ylim([-4,5])
   
        text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
            ['$\beta_{\mathrm{RA}} = ', sprintf('%3.2f',Fit.RA.beta(1)),'$'])
        text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
            ['$\beta_{\mathrm{',plot_method,'}} = ', sprintf('%3.2f',Fit.(plot_method).beta(1)),'$'])
    
       set(gca,'TickLabelInterpreter', 'LaTeX')
end




    