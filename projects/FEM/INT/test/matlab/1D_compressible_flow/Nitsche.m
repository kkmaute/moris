%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%% generate Nitsche contributions for 1D compressible flow model

close all;
clear;
clc;

%% Notes

% inputs needed for nitsche
% rho-, vel-, and temp-dofs
% normal (-1 or 1 )
% time slab ( t1 and t2 )
% element ( x1 and x2 )
% beta ( 1 or -1 for sym. or unsym. Nitsche )
% gamma ( Nitsche penalty parameter )

%% constants
R=287.058;       % specific gas constant for air
muref=1.716e-5;  % reference viscosity for air
Tref=273.0;     % reference tempereature for air

%% initial definition of symbolic variables

% define global coordinate
x=sym('x','real');

% define global time
t=sym('t','real');

% define element geometry
xA=sym('xA','real');
xB=sym('xB','real');

% define time slabs: current and previous
t1=sym('t1','real');
t2=sym('t2','real');

% define nitsche params
gam = sym( 'gam', 'real' );
beta = sym( 'beta', 'real' );

% define normal
normal = sym( 'normal', 'real' );

% define dofs of current time slab
% at t1: dof1 dof2   at t2: dof3 dof4
% at xA: dof1 dof3   at xB: dof2 dof4
rdof=sym('rdof',[4 1],'real');
vdof=sym('vdof',[4 1],'real');
tdof=sym('tdof',[4 1],'real');

% define natural coordinate
r=(2*x-(xA+xB))/(xB-xA);
s=(2*t-(t1+t2))/(t2-t1);

% define shape functions in space and time
Nr=0.5*[1-r 1+r];
Ns=0.5*[1-s 1+s];
N=reshape(Nr'*Ns,1,[]);
    
% interpolate state varialbes in space and time
% in current and previous time slab
rho=N*rdof;
vel=N*vdof;
tmp=N*tdof;

% spatial derivatives
dNdx=diff(N,x);
dveldx=diff(vel,x);

% for debugging purposes
symvec=[rdof'     vdof'     tdof'                 xA  xB   t    t1  t2  normal beta gam  ];
numvec=[1 1 1 1   0 0 0 0   Tref Tref Tref Tref   0.0 1.0  0.05 0.1 0.2 -1.0   1.0  100.0];

% set integration point to left of right node depending on normal direction
x_expr = 0.5 * ( xB*(normal+1) - xA*(normal-1) );
xr = subs(r,x,x_expr);
xNr = subs(Nr,x,x_expr);
xN = subs(N,x,x_expr);
xrho = subs(rho,x,x_expr);
xvel = subs(vel,x,x_expr);
xtmp = subs(tmp,x,x_expr);
xdNdx = subs(dNdx,x,x_expr);
xdveldx = subs(dveldx,x,x_expr);

% residuals
rhoR = simplify( int( xN'* beta*R*xtmp*xvel*normal, t, t1, t2 ) );
velR = simplify( int( xN'* (gam*xvel - muref*xdveldx*normal/3.0 + xrho*R*xtmp*normal), t, t1, t2 ) );
tmpR = simplify( int( xN'* beta*xrho*R*xvel*normal, t, t1, t2 ) );

% jacobian contributions for density residual
rhoRrho = zeros(length(xN));
rhoRvel = simplify( int( xN' * (beta*R*xtmp*normal) * xN, t, t1, t2 ) );
rhoRtmp = simplify( int( xN' * (beta*R*xvel*normal) * xN, t, t1, t2 ) );

% jacobian contributions for velocity residual
velRrho = simplify( int( xN' * (R*xtmp*normal) * xN, t, t1, t2 ) );
velRvel = simplify( int( xN' * (gam*xN - muref/3*normal*xdNdx), t, t1, t2 ) );
velRtmp = simplify( int( xN' * (xrho*R*normal) * xN, t, t1, t2 ) );

% jacobian contributions for temperature residual
tmpRrho = simplify( int( xN' * (beta*R*xvel*normal) * xN, t, t1, t2 ) );
tmpRvel = simplify( int( xN' * (beta*xrho*R*normal) * xN, t, t1, t2 ) );
tmpRtmp = zeros(length(N));

% assemble into one block
Resvec = simplify( [rhoR;velR;tmpR] );
Jacmat = simplify( [rhoRrho,rhoRvel,rhoRtmp; ...
                    velRrho,velRvel,velRtmp; ...
                    tmpRrho,tmpRvel,tmpRtmp] );

% input dof vector
dofvec=[rdof;vdof;tdof];
                
% generate matlab function
matlabFunction(Resvec,Jacmat,'File','nitsche_func.m','Optimize',false, ...
        'Vars',{dofvec,xA,xB,t1,t2,normal,beta,gam});
    
% analytically compute the jacobian and compare
%JacmatAna = jacobian( Resvec, dofvec );    
    
% debug version
% matlabFunction(Resvec,Jacmat,JacmatAna,'File','nitsche_func.m','Optimize',false, ...
%         'Vars',{dofvec,xA,xB,t1,t2,normal,beta,gam});
    
