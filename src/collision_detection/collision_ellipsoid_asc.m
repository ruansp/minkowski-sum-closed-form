%collision_ellipsoid_asc Algebraic condition stating the contact 
%                               status between two ellpsoids.
%   STATUS = collision_ellipsoid_asc(s1, s2) 
%   outputs the contact status (e.g., separated, single point contact, and 
%   overlapped) between 2 rigid ellipsoids by solving a 4th order
%   polynomial.
% 
%       ARGUMENT DESCRIPTION:
%               s1 -- Ellipsoidal model of first object
%               s1 -- Ellipsoidal model of second object
% 
%       OUTPUT DESCRIPTION:
%                 STATUS - 0 if ellipsoids are apart and 1 if overlapped.
% 
%   Example
%   -------------
%     s1 = Ellipsoid({[1.0,2.340,1.0], [0,0], 5*rand(3,1), rand(1,4), [10,10]};
%     s1 = Ellipsoid({[1.0,1.0,1.5], [0,0], 5*rand(3,1), rand(1,4), [10,10]};
%     status = collision_ellipsoid_asc(s1, s2);
% 
% See also Ellipsoid, SuperQuadrics.

% References: 
%   Wang, W., Wang, J., Kim, M.-S.
%   An algebraic condition for the separation of two ellipsoids. 
%   Computer Aided Geometric Design,
%   18(6):531�539, 2001.
%
%   Jia, X., Choi, Y.-K., Mourrain, B., Wang, W.
%   An algebraic approach to continuous collision detection for ellipsoids.
%   Computer Aided Geometric Design,
%   28:164�176, 2011.
% 
% Credits:
% Daniel Simoes Lopes
% IDMEC
% Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% danlopes (at) dem ist utl pt
% http://web.ist.utl.pt/daniel.s.lopes/
%
% July 2011 original version.
%
% Modified by Sipu Ruan, ruansp@nus.edu.sg, 2021
% 
% September 2021 modified version: compatible with SuperQuadric and
% Ellipsoid classes
%
%__________________________________________________________________________
%  Characteristic polynomial:
%  f(lambda) = det(lambda*A - Ma'*(Mb^-1)'*B*(Mb^-1)*Ma)

function status = collision_ellipsoid_asc(s1, s2)
% Ellipsoid matrices in the canonical form.
A = [diag( s1.a.^(-2) ), zeros(3,1); 0, 0, 0, -1];
B = [diag( s2.a.^(-2) ), zeros(3,1); 0, 0, 0, -1];

% Rigid body transformations.
Ma = [quat2rotm(s1.q), s1.tc; 0, 0, 0, 1];
Mb = [quat2rotm(s1.q), s1.tc; 0, 0, 0, 1];

% aij belongs to A in det(lambda*A - Ma'*(Mb^-1)'*B*(Mb^-1)*Ma).
a11 = A(1,1);
a12 = A(1,2);
a13 = A(1,3);
a14 = A(1,4);
a21 = A(2,1);
a22 = A(2,2);
a23 = A(2,3);
a24 = A(2,4);
a31 = A(3,1);
a32 = A(3,2);
a33 = A(3,3);
a34 = A(3,4);
a41 = A(4,1);
a42 = A(4,2);
a43 = A(4,3);
a44 = A(4,4);

% bij belongs to b = Ma'*(Mb^-1)'*B*(Mb^-1)*Ma .
aux = Mb \ Ma;
b = aux'*B*aux;
b11 = b(1,1);
b12 = b(1,2);
b13 = b(1,3);
b14 = b(1,4);
b21 = b(2,1);
b22 = b(2,2);
b23 = b(2,3);
b24 = b(2,4);
b31 = b(3,1);
b32 = b(3,2);
b33 = b(3,3);
b34 = b(3,4);
b41 = b(4,1);
b42 = b(4,2);
b43 = b(4,3);
b44 = b(4,4);

% Coefficients of the Characteristic Polynomial.
T4 = (-a11*a22*a33);
T3 = (a11*a22*b33 + a11*a33*b22 + a22*a33*b11 - a11*a22*a33*b44);
T2 = (a11*b23*b32 - a11*b22*b33 - a22*b11*b33 + a22*b13*b31 - ...
      a33*b11*b22 + a33*b12*b21 + a11*a22*b33*b44 - a11*a22*b34*b43 + ...
      a11*a33*b22*b44 - a11*a33*b24*b42 + a22*a33*b11*b44 - ...
      a22*a33*b14*b41);
T1 = (b11*b22*b33 - b11*b23*b32 - b12*b21*b33 + b12*b23*b31 + ...
      b13*b21*b32 - b13*b22*b31 - a11*b22*b33*b44 + a11*b22*b34*b43 + ...
      a11*b23*b32*b44 - a11*b23*b34*b42 - a11*b24*b32*b43 + ...
      a11*b24*b33*b42 - a22*b11*b33*b44 + a22*b11*b34*b43 + ...
      a22*b13*b31*b44 - a22*b13*b34*b41 - a22*b14*b31*b43 + ...
      a22*b14*b33*b41 - a33*b11*b22*b44 + a33*b11*b24*b42 + ...
      a33*b12*b21*b44 - a33*b12*b24*b41 - a33*b14*b21*b42 + ...
      a33*b14*b22*b41);
T0 = (b11*b22*b33*b44 - b11*b22*b34*b43 - b11*b23*b32*b44 + ...
      b11*b23*b34*b42 + b11*b24*b32*b43 - b11*b24*b33*b42 - ...
      b12*b21*b33*b44 + b12*b21*b34*b43 + b12*b23*b31*b44 - ...
      b12*b23*b34*b41 - b12*b24*b31*b43 + b12*b24*b33*b41 + ...
      b13*b21*b32*b44 - b13*b21*b34*b42 - b13*b22*b31*b44 + ...
      b13*b22*b34*b41 + b13*b24*b31*b42 - b13*b24*b32*b41 - ...
      b14*b21*b32*b43 + b14*b21*b33*b42 + b14*b22*b31*b43 - ...
      b14*b22*b33*b41 - b14*b23*b31*b42 + b14*b23*b32*b41);
  
%%  Roots of the characteristic_polynomial (lambda0, ... , lambda4).
characteristic_polynomial = [T4 T3 T2 T1 T0]';
r = roots(characteristic_polynomial);

% Correct numerical error of real valued polynomial zeros that are 
% accompanied by complex numbers.
% for k = 1:4
% 	if (imag(r(k)) <= 10^-3)
% 		r(k) = real(r(k));
% 	end
% end

%%  Algebraic condition for contact detection between ellipsoids.
% Find complex roots.
complex_roots = [0,0,0,0];
for k = 1:4
    complex_roots(k) = ~isreal(r(k));
end

% Find the (real) negative roots.
negative_roots_ids = find( (~complex_roots').*r < 0);

% Contact detection status.
if length(negative_roots_ids) == 2
    if  r(negative_roots_ids(1)) ~= r(negative_roots_ids(2))
%         disp('Separation Condition: quadric surfaces are separated.')
        status = 0;
        return
    elseif abs(r(negative_roots_ids(1)) - r(negative_roots_ids(2))) <= 10^-3
            % r(negative_roots_ids(1)) == r(negative_roots_ids(2))
%         disp(['Separation Condition: quadric surfaces share a single',... 
%               'contact point.'])
        status = 1;
        return
    end
else
%     disp('Separation Condition: quadric surfaces are not separated (overlapping).')
    status = 1;
    return
end