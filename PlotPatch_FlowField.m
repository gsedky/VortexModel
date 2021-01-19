%% plot patches to delineate gust region for the plotted flowfields but use gray patch
% also delinaiate the gust profile
% W: gust width
% C: chord length
% set the limits of the patch
if ii~=1
delete(findobj('type', 'patch'));
delete(p1)
end
ylimit = ylim;
xlimit = xlim;
% Gray patch
gray = [0.04 0.04 0.04];
p = patch([-X_initial+U(ii)*Time(ii)-GustWidth -X_initial+U(ii)*Time(ii) -X_initial+U(ii)*Time(ii) -X_initial+U(ii)*Time(ii)-GustWidth],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],gray);
p.FaceAlpha=0.2;
p.EdgeColor = 'none';

%% plot simple gust profile 
x_gust = [-X_initial+U(ii)*Time(ii)-GustWidth, -X_initial+U(ii)*Time(ii)]; % x boundaries of gust
x = linspace(xlimit(1), xlimit(2), 1000);  % x vector for the flowfield window
if strcmp(GustType,'Sine-square')
y_base = - 0.9*mean(abs(ylimit));  % This s the base line y position of the line for gust velocity 0
y = y_base*ones(1,length(x));   % Create the sine-squared analog
y(1,x>x_gust(1) & x<x_gust(2))=0.5*mean(abs(ylimit))*sin((x(x>x_gust(1) & x<x_gust(2))-x_gust(2))*pi/GustWidth).^2 + y_base;
p1=plot(x,y,'k','LineWidth',1);
elseif strcmp(GustType,'Top-hat')
y_base = - 0.9*mean(abs(ylimit));  % This s the base line y position of the line for gust velocity 0
y = y_base*ones(1,length(x));   % Create the top-hat analog
y(1,x>x_gust(1) & x<x_gust(2))=0.5*mean(abs(ylimit))+ y_base;
p1=plot(x,y,'k','LineWidth',1); 
end

ylimit = ylim;
