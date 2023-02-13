%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Script for Loading Geometric Jacobians for Unit Test Data

close all;
clear;
clc;

%% Input

% choose number or spatial dimensions, choose from d={2,3}
d = 3;


%% Print all Information as Required by the Unit Test

% get workspace variables from file
FileName = [ 'geom_jacs_d' num2str( d ) '.mat' ];
load( FileName );
fprintf( 'Loading parameters from: %s \n', FileName );

% check that element is not inverted
fprintf( 'det(J) = %+1.4e \n\n', det(AEval) );

% print the evaluation point
fprintf( 'The parametric evaluation point is:' );
print_matrix( XiEval, 'tXi' );

% print the evaluation point
fprintf( 'The spatial position of the evaluation point is:' );
print_matrix( XEval, 'tX' );

% print the coefficients
fprintf( 'Nodal points of the test element:' );
print_matrix( PhysNodalPoints, 'tXHat' );

% print the geometric jacobians
fprintf( 'The geometric jacobians generated are:' );
print_matrix( JaEval, 'tJa' );
print_matrix( JbEval, 'tJb' );
print_matrix( JcEval, 'tJc' );
