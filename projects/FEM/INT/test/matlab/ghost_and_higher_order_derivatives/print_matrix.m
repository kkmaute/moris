%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [] = print_matrix( aMatrix, aName )
%PRINT_MATRIX prints a matrix for copying straight to moris source code
%   aMatrix  - matrix to print
%   aName    - variable name to print it with

%% CHECK

% get the size of the matrix
MatSize = size( aMatrix );

% check that this is indeed a matrix
assert( length( MatSize ) == 2, 'print_matrix() - aMatrix is not a matrix or vector' );

% get the number of rows and columns
NumRows = MatSize( 1 );
NumCols = MatSize( 2 );

%% PRINT

fprintf( '\n%s = { \n', aName );

for i = 1:( NumRows - 1 )

    fprintf( '        { ' );
    
    for j = 1:( NumCols - 1 )
        fprintf('%+1.15e, ', aMatrix( i, j ) );
    end
    
    fprintf('%+1.15e }, \n', aMatrix( i, NumCols ) );
    
end

fprintf( '        { ' );

for j = 1:( NumCols - 1 )
    fprintf('%+1.15e, ', aMatrix( NumRows, j ) );
end

fprintf('%+1.15e } }; \n\n', aMatrix( NumRows, NumCols ) );

%% return void
end

