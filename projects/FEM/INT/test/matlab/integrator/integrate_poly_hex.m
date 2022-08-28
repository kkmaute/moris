%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

clear all;
close all;
clc;

% set number of points
tNumberOfPoints = 5;

tPath = sprintf('/home/messe/codes/MORIS/test/src/fem/data/integrate_hex_%ix%ix%i.bin', ...
tNumberOfPoints,tNumberOfPoints,tNumberOfPoints );

%%
% set order
tOrder = 2*tNumberOfPoints - 1;

% dimension
tDimension = 3;

%%
% calculate number of entries
tNumberOfCoeffs = ( tOrder+1 )^tDimension;

% reset counter
tCount = 0;

tExponents = zeros( tDimension, tNumberOfCoeffs );


% walk through Pascal's triangle
for( k = 0:tOrder )
    a = k;
    b = 0;
    
    for( j = 1:k+1)
        c = 0;
        b = j-1;
        for( i=1:j )
            tCount = tCount + 1;
            tExponents( 1, tCount ) = a;
            tExponents( 2, tCount ) = b;
            tExponents( 3, tCount ) = c;
            c = c+1;
            b = b - 1;
        end
        
        a = a - 1;
        b = b + 1;
    end
end

tN = tOrder+1;
for k = tOrder+1:(3*tOrder)
    a = k;
    b = 0;
    for( j = 1:k+1)
        c = 0;
        b = j-1;
        for( i=1:j )
            if( ( a <= tOrder && b <= tOrder && c <= tOrder ) )
                tCount = tCount + 1;
                tExponents( 1, tCount ) = a;
                tExponents( 2, tCount ) = b;
                tExponents( 3, tCount ) = c;
                
            end
            c = c+1;
            b = b - 1;
        end
        
        a = a - 1;
        b = b + 1;
    end
end

% create polynomial
syms xi eta zeta real;
p = vpa( zeros( 1, tNumberOfCoeffs ) );
for ( k = 1:tNumberOfCoeffs )
    p( k ) =   xi^tExponents( 1, k ) ...
             * eta^tExponents( 2, k ) ...
             * zeta^tExponents( 3, k );
end


%% create random coefficoents

tCoeffs = zeros( 1, tNumberOfCoeffs );

for( i = 1:tNumberOfCoeffs )
    tCoeffs( i ) = (2*rand()-1)/2;
end

% function

f = dot( tCoeffs,p );

% integral
F = vpa( int( int( int( f, xi, -1, 1), eta, -1, 1), zeta, -1, 1 ) );

%% prepare output

tData = zeros( (tDimension+1)*tNumberOfCoeffs + 3, 1 );

% first value : number of dimensions
tData( 1 ) = tDimension;

% second value : number of coeffs
tData( 2 ) = tNumberOfCoeffs;

% third value : solution
tData( 3 ) = F;

tCount = 3;

% write coefficients
for( k=1:tNumberOfCoeffs )
    tCount = tCount + 1;
    tData( tCount ) = tCoeffs( k );
end

% write exponents
for( k = 1:tNumberOfCoeffs )
    for( i = 1:tDimension )
        tCount = tCount + 1;
        tData( tCount ) = tExponents( i, k );
    end
end

%% save data file
save_matrix_to_binary( tData, tPath );


