% Convect vortices
% Enforce the kutta condition at the trailing edge and the leading edge or
% just the trailing edge depending on the optiong "Opt"
function [sX, sY, sG] = Shedding(sX, sY, sG, bX, bY, bG, su, sv, bu, bv, dt,Opt)
% Propagate the shed vortices and Shed vortex at trailing edge   
% Note: shed wake vector [last shed from LE; past wake; last shed from TE]
if strcmp('On', Opt)
    sX = [bX(1) + bu(1)*dt;sX + su*dt; bX(end) + bu(end)*dt];                                   % The x positons of all the control points
    sY = [bY(1) + bv(1)*dt; sY + sv*dt; bY(end) + bv(end)*dt];                                   % The y positons of all the control points
    sG = [bG(1);sG; bG(end)];
else
    % Propagate the shed vortices and Shed vortex at trailing edge
sX = [sX + su*dt; bX(end) + bu(end)*dt];                                   % propagate x position and find the x position of shed vortex
sY = [sY + sv*dt; bY(end) + bv(end)*dt];                                   % propagate y position and find the y position of shed vortex
sG = [sG; bG(end)];  
end
end