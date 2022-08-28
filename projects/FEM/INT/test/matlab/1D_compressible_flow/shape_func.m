%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%clear;
%clc;

% define global coordinate
x=sym('x','real');
y=sym('y','real');

% define global time
t=sym('t','real');

% define element geometry
xA=sym('xA','real');
xB=sym('xB','real');
yA=sym('yA','real');
yB=sym('yB','real');

% define time slabs: current and previous
t1=sym('t1','real');
t2=sym('t2','real');

% parent coords
q=(2*x-(xA+xB))/(xB-xA);
r=(2*y-(yA+yB))/(yB-yA);
s=(2*t-(t1+t2))/(t2-t1);

% shape func
Nq = 0.5*[1-q 1+q];
Nr = 0.5*[1-r 1+r];
Ns = 0.5*[1-s 1+s];

N = [ Nq(1)*Nr(1)*Ns(1), ...
      Nq(2)*Nr(1)*Ns(1), ...
      Nq(2)*Nr(2)*Ns(1), ...
      Nq(1)*Nr(2)*Ns(1), ...      
      Nq(1)*Nr(1)*Ns(2), ...
      Nq(2)*Nr(1)*Ns(2), ...
      Nq(2)*Nr(2)*Ns(2), ...
      Nq(1)*Nr(2)*Ns(2) ];
 
dNdx = simplify( diff(N,x) );
dNdy = simplify( diff(N,y) );
dNdt = simplify( diff(N,t) );

B = [ dNdx , dNdx ; dNdy , dNdy ];

a = 0.211325;
b = 0.788675;
symvec = [ x y t xA xB yA yB t1 t2 ];
numvec = [ a a a 0  1  0  1  0  1  ];

BT_eval = double( subs( B', symvec, numvec ) )









