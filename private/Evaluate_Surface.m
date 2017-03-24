%
% File: Evaluate_Surface.m
% Author: Oskar Skrinjar
% Date: June 2014
% -------------------------------------------------------------------
% This function evaluates a surface at given unit vector u as f(u)*u, where
% f(u) is a function of the form
%
% f(u) = a0 + a1*Psi(u.u1) + a2*Psi(u.u2) ... + aN*Psi(u.uN).
%
% The input arguments a0, a = (a1, ..., aN), and ui (a 3 x N matrix) are
% the parameters of the surface function. The surface function is evaluated
% at the directions specified by unit vectors u, which is a 3 x M matrix.
% The output s is a 3 x M matrix contining the 3D coords of the points on
% the surface that correspond to the unit vectors specified by the u
% matrix.
%
function s = Evaluate_Surface(a0, a, ui, u)

    f = Evaluate_Surface_Function(a0, a, ui, u);
    
    s = [ f.*u(1,:) ; f.*u(2,:) ; f.*u(3,:) ];