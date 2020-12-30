function V = PolygonGen1D(PatternIn,dx,dy)
% Define polygon vertexes in [x-y] plane for FDTD simulation
% dx: grid size in nm
% dy: layer height in nm
% input unit in nm
% output unit in m

SubThickness = 7000; % in nm
N = length(PatternIn);
V = zeros(length(PatternIn)+2,2);

% Substrate vertexes
V(1,:) = [0,-SubThickness];
V(2,:) = [(N-1)*dx,-SubThickness];
% Pattern vertexes
V(2+(1:N),1) = (0:1:N-1)*dx;
V(2+(1:N),2) = PatternIn*dy;
% convert nm to m
V = V*1e-9;
end