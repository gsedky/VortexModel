
function [c,vb,alpha] = SinusoidalPitching(alpha,Time,c,vb, omega_sin,a,PitchAmplitude,U)
% Pitching motion, simple rotation
omega_sin = omega_sin*(pi/180);                                                         % pitching frequency
PitchAmplitude = PitchAmplitude*pi/180;
alpha = alpha+PitchAmplitude*sin(omega_sin*Time);                                       % the alpha time history

% body frame, velocity of contact points due to pitching in the body frame
omega = omega_sin*PitchAmplitude*cos(omega_sin*Time);                                % Instantaenouse pitch rate
for ii=1:length(U)
c.Vp_b{ii} = -omega(ii)*(c.X{1}-a);
c.Up_b{ii} = zeros(length(c.Vp_b{1}),1);
% velocity of bound vortices due to pitching in the body frame
vb.Vp_b{ii} = -omega(ii)*(vb.X{1}-a);
vb.Up_b{ii} = zeros(length(vb.Vp_b{1}),1);
end

% reshape
c.Vp_b=c.Vp_b';
c.Up_b=c.Up_b';
vb.Vp_b=vb.Vp_b';
vb.Up_b=vb.Up_b';
end