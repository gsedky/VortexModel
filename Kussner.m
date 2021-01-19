function [vg_Y,vg_X,vg_G,GR] = Kussner(C,d,GustWidth,U_final,X_initial)
% Kussner Gust
Gamma_Shear = 0.482;                           % Total vorticity in one shear layer
W = GustWidth;                               % Width of the gust
H = 40*C;                                    % height of gust
% The y and X positons of gust vortices including both shear layers
vg_Y = [-H/2:d:H/2,-H/2:d:H/2]';
vg_X = [X_initial*ones(1,length(vg_Y)/2),(X_initial-W)*ones(1,length(vg_Y)/2)]';                              % The x positons gust vortices
% The magnitude of the circulation of each vortex in the gust shear
% layer
Gamma_vortex = Gamma_Shear/(length(vg_Y)/2);
vg_G =  [(-Gamma_vortex)*ones(1,length(vg_Y)/2),(Gamma_vortex)*ones(1,length(vg_Y)/2)]';                        % Strength of vortices in the gust

% Fing the vertical velocity of the gust
[~,  V_gust] = vels(((vg_X(1)+vg_X(end))/2),0,vg_X,vg_Y,vg_G);
GR = V_gust/U_final;
end
