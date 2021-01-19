function [x, y] = Trans_bi(t, n, alpha)
% Go from body frame to inertial frame
R = [cos(alpha), sin(alpha);...
    -sin(alpha), cos(alpha) ];
% Transform shed vortex positions to the body frame
x = R(1,1)*t + R(1,2)*n;
y = R(2,1)*t + R(2,2)*n;
end