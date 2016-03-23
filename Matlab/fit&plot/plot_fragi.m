% Plot figures in untransformed space
%

function plot_fragi(Fit,Options,Pos,plot_methods)


% ------------------------------------------------
% PLOT THE NORMAL FIGURE
% ------------------------------------------------

n_methods = length(plot_methods);
n_col = Options.subsize;
n_row = ceil(n_methods/n_col);

figure('Position', [400, 200, n_col*200, n_row*400])


for ii = 1:n_methods
    
    plot_method = plot_methods{ii};
    intensity   = Fit.Common.intensity;
    eintensity  = Fit.Common.eintensity;
%     scale      = Fit.Common.scale;
    p_ci       = Fit.Common.p_ci;
    cdf_ra     = Fit.Common.cdf_ra;
    fit_cdf    = Fit.(plot_method).fragi_cdf(:,1);

    
    % fragility curve
    axes('position',Pos{ii},'XGrid','off','XMinorGrid','off','FontSize',Options.fontsize,'Box','on','Layer','top');
    
        %..........................................................................
        % IMPORTANCE CURVE AND RANGE
        %..........................................................................
        if Options.imp_range == 1
             
            importance   = Fit.RA.dpf/Fit.RA.pf;

            importance_cdf = cumtrapz(intensity, importance);
            [~, ia, ~] = unique(importance_cdf);
            ci_lower = interp1(importance_cdf(ia), intensity(ia), (1-p_ci)/2);
            ci_upper = interp1(importance_cdf(ia), intensity(ia), (1+p_ci)/2);
            
            hold on
            patch([ci_lower, ci_lower, ci_upper, ci_upper], [0, 1, 1, 0], zeros(4,1),...
                                            'EdgeColor', 'none',...
                                            'FaceColor', [1, 0.8, 0.8])
            set(gca, 'Layer','top')

             
         end
    
        plot(intensity, cdf_ra, 'black');
        hold on
        plot(eintensity, fit_cdf, 'black--');

        xlim(Options.xlim)
        text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.1*diff(ylim),...
            ['$\beta_{\mathrm{RA}} = ', sprintf('%3.2f',Fit.RA.beta(1)),'$'])
        text(min(xlim) + 0.05*diff(xlim), max(ylim) - 0.2*diff(ylim),...
            ['$\beta_{\mathrm{',plot_method,'}} = ', sprintf('%3.2f',Fit.(plot_method).beta(1)),'$'])
        
        % Add bounds
        if strcmp(plot_method,'F0')==0
            [~, hb] = boundedline(eintensity, Fit.(plot_method).median, [Fit.(plot_method).median - Fit.(plot_method).ci(:,1),  Fit.(plot_method).ci(:,2) - Fit.(plot_method).median]);
            uistack(hb, 'bottom')
        end

        title(plot_method)
        if Pos{ii}(1)<=0.15
            ylabel('$P_\mathrm{f}(s)$')
        else
            set(gca,'YTickLabel','')
        end
        xlabel('Ground snow $s\,[kN/m^2]$')
        xlim(Options.xlim)

        set(gca,'TickDir','out')
        set(gca,'TickLabelInterpreter', 'LaTeX')
        
       
     
    % density curve
    if ~strcmp(plot_method,'F0')
   
        axes('position',Pos{100+ii},'XGrid','off','XMinorGrid','off','FontSize',Options.fontsize,'Box','on','Layer','top');
          
        if median(Fit.(plot_method).pf_ratio) < 1
            ss = 1./Fit.(plot_method).pf_ratio;
            tag = ['$P_\mathrm{f,RA}/P_\mathrm{f,',plot_method,'}$'];
        else
            ss = Fit.(plot_method).pf_ratio;
            tag = ['$P_\mathrm{f,',plot_method,'}/P_\mathrm{f,RA}$'];
        end
        [~, density, x] = kde(ss);
        plot(x, density, 'black')
        hold on
        area(x, density,0, 'FaceColor', 0.9*ones(1,3))
        xlabel(tag)
        
        set(gca,'YTickLabel','')
        ylabel('$density$');
%       ylab=get(gca,'YLabel');
%       set(ylab,'Position',get(ylab,'Position' ) + [0.03 0 0])
        
        mea = mean(ss);
        med = median(ss);
        
        xlim([0,10])
        ylim([0, max(density)*1.1])
        
        plot([mea, mea], ylim, 'black')
        plot([med, med], ylim, 'black')
        ht = text(mea+0.05*diff(xlim), min(ylim)+0.1*diff(ylim), ['$mean = ', sprintf('%3.2f',mea), '$']);
        ht.Rotation = 90;
        ht = text(med+0.05*diff(xlim), min(ylim)+0.1*diff(ylim), ['$median = ', sprintf('%3.2f',med), '$']);
        ht.Rotation = 90;
        

        set(gca,'TickDir','out')
        set(gca,'TickLabelInterpreter', 'LaTeX')
    end
end




    