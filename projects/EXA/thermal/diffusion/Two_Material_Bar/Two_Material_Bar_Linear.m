%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%
% 1D FE solution (body fitted, linear elements) 
% for two material heated bar example
%

L=[2.5,4.1-2.5,5-4.1,2.5,2.5];
b=[1,1,2,2,2];
c=[1,1,0.125,0.125,0.125];
edof=[1 2;2 3;3 4;4 5;5 6]

k=zeros(6,6);
rhs=zeros(6,1);

for ie=1:5
ed=edof(ie,:);
k(ed,ed)=k(ed,ed)+c(ie)/L(ie)*[1 -1;-1 1];
rhs(ed)=rhs(ed)+0.5*b(ie)*L(ie);
end
rhs(6)=rhs(6)+5;

u=k(2:6,2:6)\rhs(2:6)
