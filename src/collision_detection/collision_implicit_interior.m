function [status, dist, pt_cls] = collision_implicit_interior(s1, s2)
% collision_implicit_interior answer proximity queries between two
% superquadrics, using implicit functions and interior-point optimization
%
%  Inputs:
%    s1, s2 : SuperQuadrics class
%
%  Outputs:
%    status : status of contact
%               true -- in collision
%               false -- separated
%    dist   : Distance between the two objects
%    pt_cls : Witness points on s1 and s2
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021
%
%  See also
%    fmincon

% Optimization
option = optimoptions('fmincon', 'Algorithm', 'interior-point',...
    'display', 'none');
S = fmincon(@(x) func(x), zeros(6,1), [], [], [],...
    [], [], [], @(x) nlcon(x, s1, s2), option);

% Distance and closest points
pt_cls.s1 = S(1:3);
pt_cls.s2 = S(4:6);
dist = norm(pt_cls.s2 - pt_cls.s1);

% Collision status
if dist > 1e-3
    status = 0;
else
    status = 1;
end
end

%% Objective function
function F = func(x)
F = norm(x(1:3) - x(4:6))^2;
end

function [c, ceq] = nlcon(x, s1, s2)
x_canonical = zeros(length(x),1);
x_canonical(1:3) = quat2rotm(s1.q)' * (x(1:3)-s1.tc);
x_canonical(4:6) = quat2rotm(s2.q)' * (x(4:6)-s2.tc);

ceq = [];
c = [s1.GetImplicitFunction(x_canonical(1:3));
    s2.GetImplicitFunction(x_canonical(4:6))];
end