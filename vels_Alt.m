%% Custom Functions
% Function to find the velocity induced at arbitrary points in the domain
% due to the location and strength of vortices
% keep induced velocity below critical radius constant
%Hi
function[u,  v] = vels_Alt(x,y,xv,yv,Gamma,r_c)
% Intialize u and v
% u and v are both row vectors with each column corresponding to a location
% in the domain
u = zeros(size(x));
v = u;
for i = 1:length(xv)
    sigma = r_c;
    r2 = (x-(xv(i))).^2+(y-(yv(i))).^2;
    a = r2 < sigma^2;
    r2(a) = 10000000000000;
    % Calculate the velocity induced at all locations due to a vortex i
    ui = -Gamma(i)/(2*pi)*(y- (yv(i)))./r2;
    vi = Gamma(i)/(2*pi)*(x- (xv(i)))./r2;
    if isnan(ui)
        ui = 0;
    elseif isnan(vi)
        vi = 0;
    end
    v = v + vi;
    u = u + ui;  
end
end

