%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% Script for Generating the Geometric Jacobians Needed for 3rd Order Derivatives

close all;
clear;
clc;

%% Input

% choose number or spatial dimensions, choose from d={2,3}
d = 2;

% choose polynomial order of multi-variate function, choose from p={1,2,3}
p = 3;


%% Generate Multi-Variate Analytic Functions

% define the parametric and physical coordinates
xi = sym( 'xi', [ d, 1 ], 'real' ); % parametric coordinates
x  = sym(  'x', [ d, 1 ] ); % spatial coords as direct expresssion
X  = sym(  'X', [ d, 1 ] ); % spatial coords as function of parent coords

% generate a geometry mapping
aG = cell( d, 1 );
for i = 1:d
    aG{i} = rand( d, p + 1 );
    [ X(i), x ] = get_multi_variate_polynomial( aG{i} );
end
X = subs( X, x, xi );

% generate a field function
aU  = rand( d, p + 1 );
[ u, x ] = get_multi_variate_polynomial( aU );

% choose a random parametric point for evaluation
XiEval = 2 * rand( d, 1 ) - 1;

% get the corresponding point in physical space
XEval = eval( subs( X, xi, XiEval ) );


%% Evaluate Nodal Point Positions

% get nodal points
[ ParamNodalPoints ] = get_nodal_points( d, 3 );
NumNodalPoints = length( ParamNodalPoints( :, 1 ) ); 

% initialize array storing physical nodal points
PhysNodalPoints = zeros( NumNodalPoints, d );

% compute the position of the nodal points
for iNode = 1:NumNodalPoints
    ParamNodalPoint = ParamNodalPoints( iNode, : )';
    PhysNodalPoint = eval( subs( X, xi, ParamNodalPoint ) );
    PhysNodalPoints( iNode, : ) = PhysNodalPoint';
end

%% Compute the Derivatives 

% initialize geometric derivatives
A = sym( 'A', [ d, d ] );
B = sym( 'B', [ d, d, d ] );
C = sym( 'C', [ d, d, d, d ] );

% initialize the arrays storing the evaluated derivatives
AEval = zeros( d, d );
BEval = zeros( d, d, d );
CEval = zeros( d, d, d, d );

% construct the derivatives of the geometry mapping
for l = 1:d    
    for i = 1:d
        A(l,i) = diff( X(l), xi(i) ); 
        AEval(l,i) = eval( subs( A(l,i), xi, XiEval ) );
        for j = 1:d
            B(l,i,j) = diff( A(l,i), xi(j) );
            BEval(l,i,j) = eval( subs( B(l,i,j), xi, XiEval ) );
            for k = 1:d
                C(l,i,j,k) = diff( B(l,i,j), xi(k) );
                CEval(l,i,j,k) = eval( subs( C(l,i,j,k), xi, XiEval ) );
            end
        end
    end  
end

% % initialize function derivatives
% P = sym( 'P', [ d, 1 ] );
% Q = sym( 'Q', [ d, d ] );
% R = sym( 'R', [ d, d, d ] );

% % construct the field function derivatives
% for l = 1:d
%     P(l) = diff( u, x(l) );
%     for m = 1:d
%         Q(l,m) = diff( P(l), x(m) );
%         for n = 1:d
%             R(l,m,n) = diff( Q(l,m), x(n) );
%         end
%     end
% end

%% Compute the Symbolic Expressions for the Geometric Jacobians

% get the number of condensed indices
[ ~, NumCondensedDoubleIndices  ] = convert_to_condensed_index( [1,1], d );
[ ~, NumCondensedTrippelIndices ] = convert_to_condensed_index( [1,1,1], d );

% get the number of double and trippel indices
[ ~, NumDoubleIndices  ] = convert_to_double_index( 1, 1, d );
[ ~, NumTrippelIndices ] = convert_to_trippel_index( 1, 1, 1, d );

% build a matrix that converts a vector of condensed indices into unique permutations for Ja
CondenseOperatorForJa = zeros( NumCondensedTrippelIndices, NumTrippelIndices );
for Row = 1:NumCondensedTrippelIndices
    Indices = convert_from_condensed_index( Row, 3, d );
    Col = convert_to_trippel_index( Indices(1), Indices(2), Indices(3), d );
    CondenseOperatorForJa(Row,Col) = 1;
end

% build a matrix that converts a vector of condensed indices into unique permutations for Jb
CondenseOperatorForJb = zeros( NumCondensedDoubleIndices, NumDoubleIndices );
for Row = 1:NumCondensedDoubleIndices
    Indices = convert_from_condensed_index( Row, 2, d );
    Col = convert_to_double_index( Indices(1), Indices(2), d );
    CondenseOperatorForJb(Row,Col) = 1;
end

% initialize the geometric jacobians
Ja = sym( 'Ja', [ NumTrippelIndices, NumCondensedTrippelIndices ] );
Jb = sym( 'Jb', [ NumTrippelIndices, NumCondensedDoubleIndices ] );
Jc = sym( 'Jc', [ NumTrippelIndices, d ] );

% initialize the arrays storing the evaluated derivatives
JaEval = zeros( NumTrippelIndices, NumCondensedTrippelIndices );
JbEval = zeros( NumTrippelIndices, NumCondensedDoubleIndices );
JcEval = zeros( NumTrippelIndices, d );

% fill the first geometric jacobians
for ijk = 1:NumTrippelIndices
    
    % get [i,j,k] for the current trippel index
    [ i, j, k ] = convert_from_trippel_index( ijk, d );
    
    % fill the first geometric jacobian
    for lmn = 1:NumCondensedTrippelIndices
        
        % get [l,m,n] for the current condensed index
        Indices = convert_from_condensed_index( lmn, 3, d );
        l = Indices( 1 );
        m = Indices( 2 );
        n = Indices( 3 );
        
        % compute the jacobian entry
        Ja( ijk, lmn ) = A(l,i) * A(m,j) * A(n,k);
        JaEval( ijk, lmn ) = AEval(l,i) * AEval(m,j) * AEval(n,k);
    end
    
    % fill the second geometric jacobian
    for lm = 1:NumCondensedDoubleIndices
        
        % get [l,m] for the current condensed index
        Indices = convert_from_condensed_index( lm, 2, d );
        l = Indices( 1 );
        m = Indices( 2 );
        
        % compute the jacobian entry
        Jb( ijk, lm ) = A(l,i)*B(m,j,k) + A(m,j)*B(l,i,k) + A(m,k)*B(l,i,j);
        JbEval( ijk, lm ) = AEval(l,i)*BEval(m,j,k) + AEval(m,j)*BEval(l,i,k) + AEval(m,k)*BEval(l,i,j);
    end
    
    % fill the third geometric jacobian
    for l = 1:d 
        
        % compute the jacobian entry
        Jc( ijk, l ) = C(l,i,j,k);
        JcEval( ijk, l ) = CEval(l,i,j,k);
    end
    
end

% condense into the needed format
JaEval = JaEval * CondenseOperatorForJa;
JbEval = JbEval * CondenseOperatorForJb;

%% Print all Information as Required by the Unit Test

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

% save current workspace variables to file for reuse
FileName = [ 'geom_jacs_d' num2str( d ) '.mat' ];
save( FileName );
fprintf( 'Saving parameters to: %s \n', FileName );
