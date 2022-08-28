%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%
% Analytical solution for two material heated bar example
%

clear all

syms x

L1=sym('L1','real');
L2=sym('L2','real');
a1=sym('a1','real');
b1=sym('b1','real');
c1=sym('c1','real');
a2=sym('a2','real');
b2=sym('b2','real');
c2=sym('c2','real');
k1=sym('k1','real');
k2=sym('k2','real');
q1=sym('q1','real');
q2=sym('q2','real');
p2=sym('p2','real');

Tprec=sym('Tpre','real');

T1=a1+b1*x+c1*x^2;
T2=a2+b2*x+c2*x^2;

gradT1=diff(T1,x);
gradT2=diff(T2,x);

grad2T1=diff(gradT1,x);
grad2T2=diff(gradT2,x);

Eq(1)=+k1*grad2T1+q1;              % bulk
Eq(2)=+k2*grad2T2+q2;              % bulk
Eq(3)=subs(T1,x,0)-Tprec;          % Dbc
Eq(4)=k2*subs(gradT2,x,L1+L2) - p2;% Nbc
Eq(5)=subs(T1,x,L1)-subs(T2,x,L1); % Ifc - temperature continuity
Eq(6)=k1*subs(gradT1,x,L1)-k2*subs(gradT2,x,L1); % Ifc - flux continuity

for i=1:6
    R(i,1)=subs(Eq(i),[a1,b1,c1,a2,b2,c2],[0,0,0,0,0,0]);
    J(i,1)=diff(Eq(i),a1);
    J(i,2)=diff(Eq(i),b1);
    J(i,3)=diff(Eq(i),c1);
    J(i,4)=diff(Eq(i),a2);
    J(i,5)=diff(Eq(i),b2);
    J(i,6)=diff(Eq(i),c2);
end
sol=-simplify(J\R);

Tsol1=simplify(subs(T1,[a1,b1,c1],sol(1:3)'));
Tsol2=simplify(subs(T2,[a2,b2,c2],sol(4:6)'));

Tsol1
Tsol2
diff(Tsol1,x)
diff(Tsol2,x)
