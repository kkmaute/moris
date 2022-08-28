%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%
% 1D FE solution (body fitted, quadratic elements) 
% for two material heated bar example
%

L=[2.5,4.1-2.5,5-4.1,2.5,2.5];
b=[1,1,2,2,2];
c=[1,1,0.125,0.125,0.125];
edof=[1 2 3;3 4 5;5 6 7;7 8 9;9 10 11]

k=zeros(11,11);
rhs=zeros(11,1);

for ie=1:5
ed=edof(ie,:);
k(ed,ed)=k(ed,ed)+c(ie)/L(ie)/3*[7 -8 1;-8 16 -8;1 -8 7];
rhs(ed)=rhs(ed)+b(ie)*L(ie)/6*[1; 4; 1];
end
rhs(11)=rhs(11)+5;

u=k(2:11,2:11)\rhs(2:11)
