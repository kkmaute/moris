%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

function [ GaussPoints, Weights ] = get_gauss_points( d, nq )
%GET_GAUSS_POINTS Get the gauss points and weights for a unit cube [0,1]^d
%   d  - number of spatial dimensions
%   nq - number of quadrature points per spatial dimension

%% Generate Evaluation Points

% 1D Gauss point rule (on [-1,1] interval)
if ( nq == 1 )
    GPs1D = [ 0 ];
    Ws1D = [ 2 ];
elseif ( nq == 2 )
    GPs1D = [ -1/sqrt(3), 1/sqrt(3) ];
    Ws1D = [ 1, 1 ];
elseif ( nq == 3 )
    GPs1D = [ -sqrt(3/5), 0, sqrt(3/5) ];
    Ws1D = [ 5/9, 8/9, 5/9 ];
else
    assert( false, 'get_gauss_points() - Number of quadrature points not supported.' );
end

% transform onto [0,1] interval
GPs1D = 0.5 * GPs1D + 0.5;
Ws1D = 0.5 * Ws1D;

% get the number of GPs
NumGPs = length( GPs1D );

% initialize arrays of Gauss points and weights
GaussPoints = zeros( NumGPs^d, d );
Weights = ones( NumGPs^d, 1 );

% index which GP combination is used
IndVec = ones( 1, d );
Complete = false;
Count = 0;

while ~Complete
    
    % count the number of Gauss points that have been added
    Count = Count + 1;
    
    % fill point and compute weight for current Gauss Point
    for i = 1:d
        GaussPoints( Count, i ) = GPs1D( IndVec( i ) );
        Weights( Count ) = Weights( Count ) * Ws1D( IndVec( i ) );
    end
    
    % check if this is the last gauss point to be added
    Complete = true;
    for i = 1:d
        if( IndVec(i) < NumGPs )
            Complete = false;
        end
    end
    
    % check that the completion of the loop happens when expected
    if( Complete )
        assert( Count == NumGPs^d, 'get_gauss_points() - Gauss point generation loop failed.' );
        break;
    end
    
    % increment vector indicating which GP coordinates to use
    LoopBroke = false;
    for i = 1:d
        if( IndVec(i) < NumGPs )
            IndVec(i) = IndVec(i) + 1;
            LoopBroke = true;
            break;
        else
            IndVec(i) = 1;
        end
    end
    assert( LoopBroke, 'get_gauss_points() - Something went wrong in the Gauss generation point loop.' );

end

%% return GaussPoints, Weights
end

