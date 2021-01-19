% Shed a vortex from the trailing edge and propagate shed vortices
% Enforce the kutta condition at the trailing edge
function [sX, sY, sG] = TE_Shedding(sX, sY, sG, bX, bY, bG, su, sv, bu, bv, dt)
% Propagate the shed vortices and Shed vortex at trailing edge
sX = [sX + su*dt; bX(end) + bu(end)*dt];                                   % propagate x position and find the x position of shed vortex
sY = [sY + sv*dt; bY(end) + bv(end)*dt];                                   % propagate y position and find the y position of shed vortex
sG = [sG; bG(end)];                                                        % Strength of vortices
% Impose Kutta condition on the trailing edge of the wing
%bG(end) = 0;
end