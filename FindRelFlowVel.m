%% Find the velocity at the control points and vortices
% Find the relative flow velocity at points on the wing due 
% as well as the wing's motion, in the body frame of the wing 
% NOTE that this is the velocity of the fluid relative to the points on the
% wing
% p indicates some arbitrary point of interest on the wing, that can be a
% bound vortex or a control point
function [p_u_b, p_v_b] = FindRelFlowVel(pX,pY, vX, vY, vG,r_c,U,V,alpha,Vpb)
[p_u,  p_v] = velsLO(pX,pY, vX, vY, vG,r_c);
% add the kinematic velocity constraint due to surge and pitch
% Add the surge and plunge components in the inertial frame
p_u = p_u + U;
p_v = p_v + V;

% Switch to body frame based on alpha and add the pitching component
[p_u_b, p_v_b] = Trans_ib(p_u, p_v, alpha);           % transform
p_v_b = p_v_b - Vpb;                            % add pitching component
end