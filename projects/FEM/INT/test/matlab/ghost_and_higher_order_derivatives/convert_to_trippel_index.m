%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ ijk, n ] = convert_to_trippel_index( i, j, k, d )
%CONVERT_TO_TRIPPEL_INDEX Summary of this function goes here
%     i - index corresponding to first spatial dimension
%     j - index corresponding to second spatial dimension
%     k - index corresponding to third spatial dimension
%     d - number of spatial dimensions
%
% OUT:
%   ijk - trippel index
%     n - maximum trippel index

%% Check Inputs

assert( i <= d &&  i > 0, 'convert_to_trippel_index() - index i exceeds bounds' );
assert( j <= d &&  j > 0, 'convert_to_trippel_index() - index j exceeds bounds' );
assert( k <= d &&  k > 0, 'convert_to_trippel_index() - index k exceeds bounds' );

assert( d == 2 || d == 3, 'convert_to_trippel_index() - unknown number of spatial dimensions' );

%% Convert

% ordering conventions
if( d == 2 ) % for 2D
    IJK_Convention = [ 111, 222, 112, 122 ];
elseif( d == 3 ) % for 3D
    IJK_Convention = [ 111, 222, 333, 112, 113, 122, 223, 133, 233, 123 ];
end

% get the number of unique indices
n = length( IJK_Convention );

% create a sorted list of unique combinations
[ IJK_Convention_Sorted, IndexPosition ] = sort( IJK_Convention );

% sort the input indices and create an index
ijk_Sorted = sort( [ i, j, k ] );
ijk_Sorted = 100 * ijk_Sorted( 1 ) + 10 * ijk_Sorted( 2 ) + ijk_Sorted( 3 );

% find them in the list
PositionsInSortedList = find( IJK_Convention_Sorted == ijk_Sorted );
assert( length( PositionsInSortedList ) == 1, ...
    'convert_to_trippel_index() - Something is wrong in the index sorting convention.' );
ijk = IndexPosition( PositionsInSortedList( 1 ) );

%% return ijk, n 
end

