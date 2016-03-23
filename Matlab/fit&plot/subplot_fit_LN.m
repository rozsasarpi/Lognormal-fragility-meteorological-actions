% Initialze the figure
%

function [Pos,plot_pics] = subplot_fit_LN(fit_methods,Options)

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
sub_y = ceil(sub_pcs/sub_x);

% Geometry ratios
plotwidth = 1;            
plotheight = 1*(sub_y/sub_x);

leftmargin   = 0.15;
rightmargin  = 0.05;
bottommargin = 0.05* plotwidth/plotheight;
topmargin    = 0.05* plotwidth/plotheight;


if Options.subsize == 2
    spacex = 0.08;
else
    spacex = 0.05;
end

spacey = 0.3/sub_y;

subxsize      = (1-leftmargin-rightmargin-spacex*(sub_x-1.0))/sub_x;
subysize      = (1-topmargin-bottommargin-spacey*(sub_y-1.0))/sub_y;

% bottom_left corner points of the subfigure in x direction (from left to right)
x_points = [leftmargin, repmat(subxsize+spacex,1,sub_x-1),subxsize, rightmargin]; 

% bottom_left corner points of the subfigure in y direction (from top to bottom)
y_points = [topmargin+subysize, repmat(subysize+spacey,1,sub_y-1), bottommargin]; 

% get the cordinate
sum_x = cumsum(x_points);
sum_y = cumsum(y_points);
sum_y = 1-sum_y;

%build the position matrix
clear Pos
k=1;
for jj=1:sub_y
    for ii=1:sub_x
        Pos{k} = [sum_x(ii) sum_y(jj) subxsize subysize];
        k=k+1;
    end
end



    