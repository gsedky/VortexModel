
function [c,vb,alpha] = NoPitching(alpha,Time,c,vb,U)
alpha = linspace(alpha,alpha,length(Time)); % the alpha time history
% body frame, velocity of contact points due to pitching in the body frame
c.Vp_b{1} = zeros(size(c.X{1}));
c.Up_b{1} = zeros(size(c.X{1}));
% velocity of bound vortices due to pitching in the body frame
vb.Vp_b{1} = zeros(size(vb.X{1}));
vb.Up_b{1} = zeros(size(vb.X{1}));

% add this to the rest of the time frames
% Assign the value of the velocity vector normal to the chord at the
% intitial time to all the times and transpose to get the right cell
% shape
c.Vp_b(1:length(U))=c.Vp_b;
c.Up_b(1:length(U))=c.Up_b;
vb.Vp_b(1:length(U))=vb.Vp_b;
vb.Up_b(1:length(U))=vb.Up_b;
c.Vp_b=c.Vp_b';
c.Up_b=c.Up_b';
vb.Vp_b=vb.Vp_b';
vb.Up_b=vb.Up_b';
end
