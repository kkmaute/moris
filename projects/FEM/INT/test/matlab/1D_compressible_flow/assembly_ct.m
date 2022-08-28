%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

% debug flag
dbg = 0;

% define constants
R=287.058;       % specific gas constant for air
gam=1.4;         % ratio of specific heat
cv=0.718e3;      % specific heat for air
kappa=24.35e-3;  % conductivity of air
muref=1.716e-5;  % reference viscosity for air
Tref=273.15;     % reference tempereature for air

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

t1prev=sym('t1prev','real');
t2prev=sym('t2prev','real');

% define penalty factors for time continuity
ctrho=sym('ctrho','real');
ctvel=sym('ctvel','real');
cttmp=sym('cttmp','real');

% define dofs of current time slab
% at t1: dof1 dof2   at t2: dof3 dof4
% at xA: dof1 dof3   at xB: dof2 dof4
rdof=sym('rdof',[4 1],'real');
vdof=sym('vdof',[4 1],'real');
tdof=sym('tdof',[4 1],'real');

% define dofs of prev time slab
rprev=sym('rprev',[4 1],'real');
vprev=sym('vprev',[4 1],'real');
tprev=sym('tprev',[4 1],'real');

% define natural coordinate
r=(2*x-(xA+xB))/(xB-xA);
s=(2*t-(t1+t2))/(t2-t1);

sprev=(2*t-(t1prev+t2prev))/(t2prev-t1prev);

% define shape functions in space and time
Nr=0.5*[1-r 1+r];
Ns=0.5*[1-s 1+s];

N=reshape(Nr'*Ns,1,[]);           % shape function in current time slab

Nsprev=0.5*[1-sprev 1+sprev];
Nprev=reshape(Nr'*Nsprev,1,[]);   % shape function in previous time slab

% interpolate state varialbes in space and time
% in current and previous time slab
rho=N*rdof;
vel=N*vdof;
tmp=N*tdof;

rhoprev=Nprev*rprev;
velprev=Nprev*vprev;
tmpprev=Nprev*tprev;

% check that dofs are correclty ordered
subs(rho,t,t1);
subs(rho,t,t2);
subs(rho,x,xA);
subs(rho,x,xB);

subs(rhoprev,t,t1prev);
subs(rhoprev,t,t2prev);

mu=muref;

% for debugging purposes
if dbg
    symvec=[rdof'     vdof'     tdof'                 x        xA  xB   t        t1  t2   QA  QB ];
    numvec=[1 1 1 1   0 0 0 0   Tref Tref Tref Tref   0.211325 0.0 1.0  0.211325 0.0 1.0  0.0 0.0];
end

% time continuity penalizing jumps in rho,vel, and tmp
% note: not clear whehter instead of temp energy should be used
Nt1=subs(N,t,t1);

Rrhotime=simplify(ctrho*int(Nt1'*(subs(rho,t,t1)-subs(rhoprev,t,t2prev)),x,xA,xB));
Rveltime=simplify(ctvel*int(Nt1'*(subs(vel,t,t1)-subs(velprev,t,t2prev)),x,xA,xB));
Rtmptime=simplify(cttmp*int(Nt1'*(subs(tmp,t,t1)-subs(tmpprev,t,t2prev)),x,xA,xB));

% vector of residuals
Resvec=simplify([Rrhotime;Rveltime;Rtmptime]);
dofvec=[rdof;vdof;tdof];
dofvecprev=[rprev;vprev;tprev];

Jacmat=simplify(jacobian(Resvec,dofvec));

% generate matlab function
matlabFunction(Resvec,Jacmat,'File','ondim_tc_func.m','Optimize',false, ...
    'Vars',{dofvec,dofvecprev,xA,xB,t1,t2,t1prev,t2prev,ctrho,ctvel,cttmp});

% debug
% matlabFunction(Resvec,Jacmat,Rrhobulk,Rvelbulk,Rtmpbulk,'File','ondim_res_dbg_func.m','Optimize',false, ...
%     'Vars',{dofvec,dofvecprev,xA,xB,t1,t2,t1prev,t2prev,ctrho,ctvel,cttmp});

% stop here
return;
