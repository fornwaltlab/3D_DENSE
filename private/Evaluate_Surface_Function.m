%
% File: Evaluate_Surface_Function.m
% Author: Oskar Skrinjar
% Date: June 2014
% -------------------------------------------------------------------
% This function evaluates the surface function f(u) at given unit vector u,
% where f(u) is a function of the form
%
% f(u) = a0 + a1*Psi(u.u1) + a2*Psi(u.u2) ... + aN*Psi(u.uN).
%
% The input arguments a0, a = (a1, ..., aN), and ui (a 3 x N matrix) are
% the parameters of the surface function. f(u) is evaluated at the
% directions specified by unit vectors u, which is a 3 x M matrix. The
% output f is a 1 x M array contining the function values that correspond
% to the unit vectors specified by the u matrix.
%
function f = Evaluate_Surface_Function(a0, a, ui, u)

    if isempty(a)
        f = a0 * ones(1,size(u,2));
    else
        f = a0 + sum(repmat(a, 1, size(u,2)).*Psi(ui'*u));
    end