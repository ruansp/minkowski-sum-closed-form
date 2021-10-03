function err = kiss_point_measure(s1, s2, M1, M2, x1, mink)
% kiss_point_measure Error metric for kissing point
%  Compute the candidate kissing point via the proposed closed-form 
%  Minkowski sums. Verify the implicit function and gradient values.
%
%  Inputs:
%    s1, s2: Geometry classes
%    M1, M2: Linear tranformation matrices
%    x1    : A point on s1 boundary
%    mink  : A set of points on Minkowski sums boundary
%
%  Output:
%    err   : Error metrics for the kissing point
%            (1) implicit function value of the candidate kissing point
%            (2) gradient error of the candidate kissing point
%
%  Author:
%    Sipu Ruan, ruansp@nus.edu.sg, 2021

% Remove NaN entries
for i = 1:size(mink,1)
    idx_nan = isnan(mink(i,:));
    mink(:,idx_nan) = [];
    x1(:,idx_nan) = [];
end

% Implicit function value and gradient of computed kissing point
x2_local = M2 \ (x1 - mink);
val_implicit_2 = s2.GetImplicitFunction(x2_local);
    
m_1 = M1' \ s1.GetGradientsFromParamTapered(M1\x1, 'cartesian');
m_2 = M2' \ s2.GetGradientsFromParamTapered(x2_local, 'cartesian');
grad_dot = diag(m_1' * m_2)';
grad_norm_prod = sqrt( sum(m_1.^2, 1) ) .* sqrt( sum(m_2.^2, 1) );

implicit_err_vec = abs(val_implicit_2);
grad_err_vec = abs(grad_dot./grad_norm_prod - cos(pi));

err(1) = mean(implicit_err_vec, 'omitnan');
err(2) = mean(grad_err_vec, 'omitnan');
end
