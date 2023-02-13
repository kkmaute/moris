%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ i, j ] = convert_from_double_index( ij, d )
%CONVERT_FROM_DOUBLE_INDEX Convert double index back to individual indices
%   ij - double index to unpack
%    d - number of spatial dimensions

%% Check Input

if( d == 2 )
    assert( ij <= 3 &&  ij > 0, 'convert_from_double_index() - trippel index exceeds bounds' );
elseif( d == 3 )
    assert( ij <= 6 &&  ij > 0, 'convert_from_double_index() - trippel index exceeds bounds' );
else
    assert( false, 'convert_from_double_index() - unknown number of spatial dimensions' );
end
    
%% Convert

% 2D
if( d == 2 )
    
    I = [ 1, 2, 1 ];
    J = [ 1, 2, 2 ];
    
    i = I(ij);
    j = J(ij);
    
    return;
    
% 3D
elseif( d == 3 )
    
    I = [ 1, 2, 3, 2, 1, 1 ];
    J = [ 1, 2, 3, 3, 3, 2 ];
    
    i = I(ij);
    j = J(ij);
    
    return;

end

%% return i, j
end

