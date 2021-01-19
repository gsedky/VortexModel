% Define function for step repsonse
function [U,V] = StepResponse(U_final,Time)
% Surge and plunging velocity velocity
U = repmat(U_final,1,length(Time));                % Linear velocity profile of the wing
V = zeros(1,length(U));
end