%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ i, j, k ] = convert_from_trippel_index( ijk, d )
%CONVERT_FROM_TRIPPEL_INDEX Convert trippel index back to individual indices
%   ijk - trippel index to unpack
%     d - number of spatial dimensions

%% Check Input

if( d == 2 )
    assert( ijk <= 4 &&  ijk > 0, 'convert_from_trippel_index() - trippel index exceeds bounds' );
elseif( d == 3 )
    assert( ijk <= 10 &&  ijk > 0, 'convert_from_trippel_index() - trippel index exceeds bounds' );
else
    assert( false, 'convert_from_trippel_index() - unknown number of spatial dimensions' );
end
    
%% Convert

% 2D
if( d == 2 )
    
    I = [ 1, 2, 1, 1 ];
    J = [ 1, 2, 1, 2 ];
    K = [ 1, 2, 2, 2 ];
    
    i = I(ijk);
    j = J(ijk);
    k = K(ijk);
    
    return;
    
% 3D
elseif( d == 3 )
    
    I = [ 1, 2, 3, 1, 1, 1, 2, 1, 2, 1 ];
    J = [ 1, 2, 3, 1, 1, 2, 2, 3, 3, 2 ];
    K = [ 1, 2, 3, 2, 3, 2, 3, 3, 3, 3 ];
    
    i = I(ijk);
    j = J(ijk);
    k = K(ijk);
    
    return;

end

%% return i, j, k
end

