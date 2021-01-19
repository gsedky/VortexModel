% Function to find the velocity induced at arbitrary points in the domain
% due to the location and strength of vortices using Lamb-Oseen model
% function[u,  v] = vels(x,y,xv,yv,Gamma)
% Input: 
% x: vector or grid of x positions of points of interest
% y: vector or grid of y positions of points of interest
% xv: vector of x positions of vortices
% yv: vector of y positions of vortices
% Gamma: vector of the vortex strengths of each vortex point
% Output:
% u: horizental velocoity at each point of interest
% v: vertical velocoity at each point of interest
function[u,  v] = velsLO(x,y,xv,yv,Gamma,r_c)
% Intialize u and v
% u and v are both row vectors with each column corresponding to a location
% in the domain
u = zeros(size(x));
v = u;

% Take care of the fact that their might be no shed vortices
for i = 1:length(xv)
    r2 = (x-(xv(i))).^2+(y-(yv(i))).^2;
    
    % Calculate the velocity induced at all locations due to a vortex i
    % Lamb-Oseen Vortex
    ui = -(Gamma(i)/(2*pi))*((y-(yv(i)))./r2).*(1-exp((-r2)./(r_c.^2)));
    vi = (Gamma(i)/(2*pi))*((x-(xv(i)))./r2).*(1-exp((-r2)./(r_c.^2)));
    
    % If the two vortices lie on top of each other, r2 is really 0.
    % Theoretically that should not be a problem for Lamb-Oseen because 
    % (1-exp(-r^2/r_c^2)) goes to 0 faster than 1/(2pi*r) but numerically
    % that expression still blows up..., NaN*0=NaN not 0
    % Correct the problem by removing the nans
    ui(isnan(ui))=0;
    vi(isnan(vi))=0;   
    v = v + vi;
    u = u + ui;
end
end