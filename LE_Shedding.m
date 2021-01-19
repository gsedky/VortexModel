% Shed a vortex from the trailing edge and the leading edge and propagate 
% shed vortices
% Enforce the kutta condition at the trailing edge and the leading edge
function [sX, sY, sG] = LE_Shedding(sX, sY, sG, bX, bY, bG, su, sv, bu, bv, dt)
% Propagate the shed vortices and Shed vortex at trailing edge   
% Note: shed wake vector [last shed from LE; past wake; last shed from TE]
    sX = [bX(1) + bu(1)*dt;sX + su*dt; bX(end) + bu(end)*dt];                                   % The x positons of all the control points
    sY = [bY(1) + bv(1)*dt; sY + sv*dt; bY(end) + bv(end)*dt];                                   % The y positons of all the control points
    sG = [bG(1);sG; bG(end)];

end