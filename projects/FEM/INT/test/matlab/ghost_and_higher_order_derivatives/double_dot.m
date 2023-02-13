%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ DoubleDotProduct ] = double_dot( A, B )
%DOUBLE_DOT Double dot product of two tensors A_ij * B_ij
%   A,B - tensors of equal dimensions

%% FUNCTION BODY

% check size
assert( sum( size( A ) == size( B ) ) == 2, ...
    'double_dot() - Dimensions of A and B do not match.' );

% get size of the arrays
ArraySize = size( A );
NumRows = ArraySize( 1 );
NumCols = ArraySize( 2 );

% perform double dot product
DoubleDotProduct = 0;
for i = 1:NumRows
    for j = 1:NumCols
        DoubleDotProduct = DoubleDotProduct + A(i,j) * B(i,j);
    end
end

%% return DoubleDotProduct
end

