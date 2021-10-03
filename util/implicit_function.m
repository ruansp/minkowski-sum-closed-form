function val = implicit_function(s1, s2, M1, M2, mink)
% implicit_function Implicit function value of the point into a
% superquadric model
%
%  Inputs:
%    s1, s2: Geometry classes
%    M1, M2: Linear tranformation matrices
%    mink  : A set of points on Minkowski sums boundary
%
%  Output:
%    val   : Implicit function value of the point set
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021

% Discrete points on the surface of S_2
x2_in_s1 = M1 \ (M2 * s2.GetPoints() + mink);

% Skip when entries in mink is NaN
is_nan = false;
for i = 1:size(mink,1)
    if isnan(mink(i))
        is_nan = true;
        break;
    end
end

if is_nan
    val = zeros(1,size(x2_in_s1,2));
    return;
end

% Distance indicator
val = s1.GetImplicitFunction(x2_in_s1);

end