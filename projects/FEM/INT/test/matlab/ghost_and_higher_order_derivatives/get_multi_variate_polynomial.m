%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ Function, x ] = get_multi_variate_polynomial( aCoefficients )
%GET_MULTI_VARIATE_POLYNOMIAL Generates a symbolic expression of a
%multi-variate polynomial of order p in d variables.
%   aCoefficients - array of coefficients of size d x p+1

%% FUNCTION BODY

% get the size of the array
SizeCoeffs = size( aCoefficients );

% get the number of spatial dimensions
d = SizeCoeffs( 1 );

% get the polynomial order
p = SizeCoeffs( 2 ) - 1;

% initialize variables
x = sym( 'x', [ d, 1 ], 'real' );

% get the single-variate components of the multi-variate functions to be constructed
f = sym( 'f', [ d, 1 ] );
for i = 1:d
    f(i) = 0;
    for n = 1:p+1
        f(i) = f(i) + aCoefficients(i,n) * x(i)^(n-1);
    end
end

% assemble the multi-variate field functions 
Function = 1;
for i = 1:d
    Function = Function * f( i );
end

%% return Function, x
end

