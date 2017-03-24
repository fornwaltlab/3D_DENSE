%
% File: Fit_Surface.m
% Author: Oskar Skrinjar
% Date: June 2014
% -------------------------------------------------------------------
% This function fits a surface of the form S = f(u)*u, where S is the
% surface vector and u is a unit vector (the model is centered at the
% coordinate origin), and f(u) is the following scalar function
%
% f(u) = a0 + a1*Psi(u.u1) + a2*Psi(u.u2) ... + aN*Psi(u.uN).
%
% The input argument u is a 3 x N matrix containing the N unit vectors. The
% input argument v is a 3 x M matrix containing points that need to be
% fitted. The input
% argument lambda is the smoothness lambda. The output argument a0 is the a0 from the
% above formula and the output argument a is an array with elements a1, a2,
% ..., aN.
%
% If the surface cannot be fitted (i.e. if the system is not full rank)
% then a0 = [] is returned.
%
function [a0, a] = Fit_Surface(u, v, lambda)

    % M, N
    M = size(v,2);
    N = size(u,2);
    
    % compute v magnitudes
    vm = sqrt(sum(v.*v));
    
    % compute v unit vectors
    ind = find(vm <= eps);
    vm(ind) = 1; % to avoid division by zero
    vu = [ v(1,:)./vm ; v(2,:)./vm ; v(3,:)./vm ];
    vm(ind) = 0;
    
    % special case (N = 0)
    if (N == 0)
        a0 = sum(vm)/M;
        a = [];
        return;
    end
    
    
    %
    % Construct the system matrix
    %
    
    % allocate space for the system matrix
    A = zeros(N+1,N+1);
    
    % element (1,1)
    A(1,1) = 1;
    
    
    % first row and column
    A(1,2:(N+1)) = sum(Psi(vu'*u))/M;
    A(2:(N+1),1) = A(1,2:(N+1))'; % due to the system matrix symmetry
    
    % the central part from 2 to N+1 in both directions
    A(2:(N+1),2:(N+1)) = Psi(u'*vu)*Psi(vu'*u)/M;
    
    % the main diagonal (except the element 1,1) has (M/N)*lambda
    ind = sub2ind(size(A), 2:(N+1), 2:(N+1));
    A(ind) = A(ind) + lambda/N;

    
    %
    % Construct the system right hand side
    %
    
    % allocate space for the system right hand side
    b = zeros(N+1,1);
    
    % element 1
    b(1) = sum(vm)/M;

    % elements 2 to N+1
    b(2:(N+1)) = sum(repmat(vm,N,1).*Psi(u'*vu), 2)/M;
    
    
    %
    % Solve the sytem
    %
    
    % check if the system is regular
    if (rank(A) < (N+1))
        fprintf('Fit_Surface: the system is not full rank! Rank(A) = %d   N+1 = %d\n', rank(A), N+1);
        a0 = [];
        a = [];
        return;
    end
    
    % solve the system
    x = A\b;
    
    
    %
    % Set the output arguments
    %
    a0 = x(1);
    a = x(2:(N+1));