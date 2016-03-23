% Initialze the figure
%

function [Pos, plot_pics] = subplot_fit(fit_methods, Options)

% ------------------------------------------------
%INITIALIZING THE FIGURE
% ------------------------------------------------

%variable:
sub_pcs     = length(fit_methods);
plot_pics   = fit_methods;
% number of subplot
if Options.F0 == 0
    sub_pcs = sub_pcs-1;
    plot_pics{1} = [];
end
sub_x = Options.subsize;
sub_y = ceil(sub_pcs/sub_x)*2;



% Geometry ratios
plotwidth   = 1;
plotheight  = 0.85*(sub_y/sub_x);

leftmargin   = 0.1;
rightmargin  = 0.05;
bottommargin = 0.1* plotwidth/plotheight;
topmargin    = 0.05* plotwidth/plotheight;

rat_dens_frag = 0.7; %ratio between the density and fragility curve height


if Options.subsize==2
    spacex = 0.08;
else
    spacex = 0.05;
end

spacey_tot      = 0.2;
rat_space_y = 0.6; % y space ratio of the gap between the (dens and frag) and farg 

spacey_new   = spacey_tot/(sub_y/2*rat_space_y+sub_y/2-1);
spacey_dens  = spacey_new*rat_space_y;

subxsize      = (1-leftmargin-rightmargin-spacex*(sub_x-1.0))/sub_x;
subysize_frag = (1-topmargin-bottommargin-spacey_tot)/sub_y*(2/(1+rat_dens_frag));
subysize_dens = (1-topmargin-bottommargin-spacey_tot)/sub_y*(2-2/(1+rat_dens_frag));

% bottom_left corner points of the subfigure in x direction (from left to right)
x_points = [leftmargin, repmat(subxsize+spacex,1,sub_x-1),subxsize, rightmargin];

% bottom_left corner points of the subfigure in y direction (from top to bottom)
y_points = [topmargin+subysize_frag, subysize_dens+spacey_dens, repmat([subysize_frag+spacey_new,subysize_dens+spacey_dens],1,(sub_y-2)/2) bottommargin];

% get the cordinate
sum_x   = cumsum(x_points);
sum_y   = cumsum(y_points);
sum_y   = 1-sum_y;

%build the position matrix
clear Pos % why?
k = 1;
kk = 101;
for jj=1:sub_y
    for ii=1:sub_x
        if mod(jj,2)==0
            Pos{kk} = [sum_x(ii) sum_y(jj) subxsize subysize_dens];
            kk=kk+1;
        else
            Pos{k} = [sum_x(ii) sum_y(jj) subxsize subysize_frag];
            k=k+1;
        end
    end
end


% f=figure('visible','on', 'PaperUnits', 'centimeters','Position',[0 0 plotwidth_px plotheight*plotwidth_px]);
% clf(f);
end


