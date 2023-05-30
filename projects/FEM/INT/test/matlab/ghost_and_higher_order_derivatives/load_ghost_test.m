%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Script for Loading Functions and Values from Ghost IWG Unit Test Data

close all;
clear;
clc;

%% Input

% choose number or spatial dimensions, choose from d={2,3}
d = 3;

% choose polynomial order of multi-variate function, choose from p={1,2,3}
p = 3;

%% Print all Information as Required by the Unit Test

% save current workspace variables to file for reuse
FileName = ['ghost_data_d' num2str( d ) '_p' num2str( p ) '.mat'];
load( FileName );

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

