% Define function for trapezoidal motion
function [U,V] = LinearRamp(U_final,Time,Acceleration)
% Uniform motion
U = U_final*ones(1,length(Time));
% Linear increase in speed in the beggining
U(Time<U_final/Acceleration) = Acceleration*Time(Time<U_final/Acceleration);
% No plunging
V = zeros(1,length(Time));

end