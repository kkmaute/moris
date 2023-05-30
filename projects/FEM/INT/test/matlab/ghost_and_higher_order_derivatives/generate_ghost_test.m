%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Script for Generating Functions and Values for the Ghost IWG Unit Test

close all;
clear;
clc;

%% Input

% choose number or spatial dimensions, choose from d={2,3}
d = 3;

% choose polynomial order of multi-variate function, choose from p={1,2,3}
p = 3;

% number of Gauss points per spatial dimensions, choose from nq={1,2,3}
nq = 2;

%% Generate the Analytic Function and its Derivatives

% generate random arrays of coefficients
aM = rand( d, p + 1 ); % leader element
aS = rand( d, p + 1 ); % follower element
% aM = ones( d, p + 1 );
% aS = 0.5 * ones( d, p + 1 );

% generate symbolic multi-variate polynomials from coefficients
[ uM, ~ ] = get_multi_variate_polynomial( aM );
[ uS, x ] = get_multi_variate_polynomial( aS );

% get the derivatives
GradUM = simplify( gradient( uM, x ) );
GradUS = simplify( gradient( uS, x ) );

% get the 2nd derivatives
if( p > 1 )
    Grad2UM = simplify( jacobian( GradUM, x ) );
    Grad2US = simplify( jacobian( GradUS, x ) );
end

% get the 3rd derivatives
if( p > 2 )
    Grad3UM = sym( 'Grad3UM', [ d, d, d ] );
    Grad3US = sym( 'Grad3US', [ d, d, d ] );
    for i = 1:d
        for j = 1:d
            Grad3UM(i,j,:) = gradient( Grad2UM(i,j), x );
            Grad3US(i,j,:) = gradient( Grad2US(i,j), x );
        end
    end
end
        

%% Generate Coefficients for field

% get the nodal points for the element
NodalPoints = get_nodal_points( d, p );
NodalPoints = 0.5 * NodalPoints + 0.5;

% get the number of nodal points
NumNodes = length( NodalPoints(:,1) );

% compute coefficients for leader and follower elements
UHatM = zeros( NumNodes, 1 );
UHatS = zeros( NumNodes, 1 );
for i = 1:NumNodes
    xEval = NodalPoints( i, : )';
    UHatM( i ) = eval( subs( uM, x, xEval ) );
    UHatS( i ) = eval( subs( uS, x, xEval ) );
end

%% Evaluation Points

% generate a normal
Normal = rand( d, 1 );
Normal = Normal / norm( Normal );

% generate evaluation points and weights
[ GaussPoints, Weights ] = get_gauss_points( d, nq );
NumGPs = length( Weights );

%% Evaluate Ghost Integrand

% initialize arrays of integrand values for all Gauss points
IntegrandComponents = zeros( NumGPs, p );
Residual = 0;

% loop over Gauss points
for iGP = 1:NumGPs
    
    % get the current Gauss point and weight
    GaussPoint = GaussPoints( iGP, : )';
    Weight = Weights( iGP );

    % evaluate jump in first derivatives and add its contribution
    Jump = Normal' * ( eval( subs( GradUM, x, GaussPoint ) ) - eval( subs( GradUS, x, GaussPoint ) ) );
    IntegrandComponents( iGP, 1 ) = Jump^2;
    Residual = Residual + Weight * IntegrandComponents( iGP, 1 );
    
    % evaluate jump in second derivatives and add its contribution
    if( p > 1 )
        Jump = Normal' * ( eval( subs( Grad2UM, x, GaussPoint ) ) - eval( subs( Grad2US, x, GaussPoint ) ) );
        IntegrandComponents( iGP, 2 ) = Jump * Jump';
        Residual = Residual + Weight * IntegrandComponents( iGP, 2 );
    end
    
    % evaluate jump in third derivatives and add its contribution
    if( p > 2 )
        NormalJump = zeros( d, d );
        for j = 1:d
            for k = 1:d
                for i = 1:d
                    Jump = eval( subs( Grad3UM( i, j, k ), x, GaussPoint ) ) - eval( subs( Grad3US( i, j, k ), x, GaussPoint ) );
                    NormalJump(j,k) = NormalJump(j,k) + Normal( i ) * Jump;
                end
            end
        end
        IntegrandComponents( iGP, 3 ) = double_dot( NormalJump, NormalJump );
        Residual = Residual + Weight * IntegrandComponents( iGP, 3 );
    end

end

% compute the Integrands for each Gauss point
Integrands = zeros( NumGPs, 1 );
for iGP = 1:NumGPs
    Integrands( iGP ) = Weights( iGP ) * sum( IntegrandComponents( iGP, : ) );
end


%% Print all Information as Required by the Unit Test

% print the normal
fprintf('\nThe generated normal is: \n');
print_matrix( Normal, 'tNormal' );

% print the coefficients
fprintf( '\nCoefficient vector for the leader element: \n' );
print_matrix( UHatM, 'u_hat_M' );
fprintf( '\nCoefficient vector for the follower element: \n' );
print_matrix( UHatS, 'u_hat_S' );

% print the result
fprintf( '\nResidual value: %.15e \n', Residual );
fprintf( '\nIntegrand values for each Gauss point: \n' );
print_matrix( Integrands, 'Integrands' );

% save current workspace variables to file for reuse
FileName = ['ghost_data_d' num2str( d ) '_p' num2str( p ) '.mat'];
save( FileName );
fprintf( '\nSaving parameters to: %s \n', FileName );
