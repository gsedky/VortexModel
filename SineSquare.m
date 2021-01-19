
function [vg_Y,vg_X,vg_G,GR,W] = SineSquare(U, X_initial)
load('GustData.mat');
load('GustVelocity.mat','v');
GR = v(length(v)/2,length(v)/2)/max(U);
% Create the x and y position grid of the vortex sheet
[Xv, Yv] = meshgrid(x,y);
vg_Y = Yv(:) - (0.5)*mean(diff(y));
vg_X = Xv(:)+X_initial;                         % Intial X position of the gust beginning from the leading edge ;
Gamma = repmat(Gamma',length(y),1);
vg_G = Gamma(:);
end