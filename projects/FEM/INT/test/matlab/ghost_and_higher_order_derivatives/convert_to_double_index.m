%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ ij, n ] = convert_to_double_index( i, j, d )
%CONVERT_TO_DOUBLE_INDEX Summary of this function goes here
%    i - index corresponding to first spatial dimension
%    j - index corresponding to second spatial dimension
%    d - number of spatial dimensions
%
% OUT:
%   ij - double index
%    n - maximum double index

%% Check Inputs

assert( i <= d &&  i > 0, 'convert_to_double_index() - index i exceeds bounds' );
assert( j <= d &&  j > 0, 'convert_to_double_index() - index j exceeds bounds' );

assert( d == 2 || d == 3, 'convert_to_double_index() - unknown number of spatial dimensions' );

%% Convert

% ordering conventions
if( d == 2 ) % for 2D
    IJ_Convention = [ 11, 22, 12 ];
elseif( d == 3 ) % for 3D
    IJ_Convention = [ 11, 22, 33, 23, 13, 12 ];
end

% get the number of unique indices
n = length( IJ_Convention );

% create a sorted list of unique combinations
[ IJ_Convention_Sorted, IndexPosition ] = sort( IJ_Convention );

% sort the input indices and create an index
ij_Sorted = sort( [ i, j ] );
ij_Sorted = 10 * ij_Sorted( 1 ) + ij_Sorted( 2 );

% find them in the list
PositionsInSortedList = find( IJ_Convention_Sorted == ij_Sorted );
assert( length( PositionsInSortedList ) == 1, ...
    'convert_to_double_index() - Something is wrong in the index sorting convention.' );
ij = IndexPosition( PositionsInSortedList( 1 ) );

%% return ijk, n 
end

