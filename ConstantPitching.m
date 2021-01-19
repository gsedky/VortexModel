
function [c,vb,alpha] = ConstantPitching(alpha,Time,c,vb,omega,a,Duration,U)
% Pitching motion, simple rotation
omega = omega*(pi/180);                           % pitch rate, rad/s
alpha = linspace(alpha,alpha+omega*Duration,length(Time)); % the alpha time history
% body frame, velocity of contact points due to pitching in the body frame
c.Vp_b{1} = -omega*(c.X{1}-a);
c.Up_b{1} = zeros(length(c.Vp_b{1}),1);
% velocity of bound vortices due to pitching in the body frame
vb.Vp_b{1} = -omega*(vb.X{1}-a);
vb.Up_b{1} = zeros(length(vb.Vp_b{1}),1);

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