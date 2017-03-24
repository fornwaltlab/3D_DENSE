%
% File: Psi.m
% Author: Oskar Skrinjar
% Date: December 2013
% -------------------------------------------------------------------
% The Psi function used in the pseudo-thin-plate-splines on the sphere.
%

function y = Psi(t)

    m = 2;

    index = find( t<1 );

    q2 = (1/m)*ones( size(t,1) , size(t,2) );

    q2(index) = (  Psi_A(t(index)).*(   12*Psi_W(t(index)).^2 - 4*Psi_W(t(index))  ) - 6 * Psi_C(t(index)).*Psi_W(t(index)) + 6*Psi_W(t(index)) + 1 )/ 2;

    y = ( 1/( factorial( 2*m-2 ) ) * q2 - 1/( factorial( 2*m-1 ) ) ) / (2*pi);


function y = Psi_A( z )

    y = log( 1 + 1./sqrt( Psi_W(z)));
    
    
function y = Psi_C( z )

    y = 2 * sqrt(  Psi_W(z) );
    
    
function y = Psi_W( z )

    y = (1 - z)/2;