function y = rescale(x,type)
% Release:
%   Q = rescale(u_norm,"release")
%
% Power:
%   P = rescale([u_norm V_norm a b],"power")

% Flow scaling parameters
std_Q = 686;
std_q = 0.004;
alpha_q = std_Q/std_q;
dt = 3600; % time step of one hour

% Head scaling parameters
H0  = 10;
dH  = 10;
a   = 5;
b   = 0.45;

% Power scaling parameters
rho = 1000;
g   = 9.81;
eta = 0.9;
F = 750;

% Hydropower conversion constant
c = eta*rho*g/1e6;

switch lower(type)

    case "release"  % m^3/s
        y = alpha_q .* x;         

    case "power"    % MWh
        y = min(c .* x, F); % apply capacity rating 

    case "head"     % m
        x_floor = max(x, 0);  % clip negative volumes to zero
        h = a .* (x_floor .^ b);
        h_norm = min(max(h ./ a, 0), 1);
        y = H0 + dH .* (h_norm - 0.5);

    case "volume"   % m3
        y = x .* alpha_q .* dt; 

    otherwise
        error("Unknown rescale type")

end

end