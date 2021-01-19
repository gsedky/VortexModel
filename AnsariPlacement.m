% Modify the positions of the nasecent leading and trailing/only trailing
% edge vortices using the Ansari way
function [A,vbX,vbY]=AnsariPlacement(A,vbX,vbY,cX,cY,vsX,vsY,r_c,alpha,Opt)
if strcmp('On', Opt)
    % Modify the position of nascent vortices and trailing and leading edge
    % leading edge
    vbX(1)= cX(1)+(1/3)*(vsX(1)-cX(1));
    vbY(1)= cY(1)+(1/3)*(vsY(1)-cY(1));
    % trailing edge
    vbX(end)= cX(end)+(1/3)*(vsX(end)-cX(end));
    vbY(end)= cY(end)+(1/3)*(vsY(end)-cY(end));
else
    % Modify the position of nascent vortices at trailing edge
    % trailing edge
    vbX(end)= cX(end)+(1/3)*(vsX(end)-cX(end));
    vbY(end)= cY(end)+(1/3)*(vsY(end)-cY(end));
end

% update the matrix A coefficients to reflect these chagnges
% Calculate the velocity induced at all locations due to a vortex i
% Lamb-Oseen Vortex
% LEV
r2 = (vbX(2:end)-(vbX(1))).^2+(vbY(2:end)-(vbY(1))).^2;
ui = -(1/(2*pi))*((vbY(2:end)-(vbY(1)))./r2).*(1-exp((-r2)./(r_c.^2)));
vi = (1/(2*pi))*((vbX(2:end)-(vbX(1)))./r2).*(1-exp((-r2)./(r_c.^2)));
A(1:end-1,1)=ui*sin(alpha)+vi*cos(alpha);
% TEV
r2 = (vbX(1:end-1)-(vbX(end))).^2+(vbY(1:end-1)-(vbY(end))).^2;
ui = -(1/(2*pi))*((vbY(1:end-1)-(vbY(end)))./r2).*(1-exp((-r2)./(r_c.^2)));
vi = (1/(2*pi))*((vbX(1:end-1)-(vbX(end)))./r2).*(1-exp((-r2)./(r_c.^2)));
A(1:end-1,end)=ui*sin(alpha)+vi*cos(alpha);
end