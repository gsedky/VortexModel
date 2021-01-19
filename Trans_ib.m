function [t, n] = Trans_ib(x, y, alpha)
% Go from inertial frame to body frame
R = [cos(alpha), -sin(alpha);...
    sin(alpha), cos(alpha) ];
% Transform shed vortex positions to the inertial frame
t = R(1,1)*x + R(1,2)*y;
n = R(2,1)*x + R(2,2)*y;
end
