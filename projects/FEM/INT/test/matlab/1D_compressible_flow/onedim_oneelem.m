%
% Copyright (c) 2022 University of Colorado
% Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%------------------------------------------------------------------------------------
%

%
% 1-D compressible flow model
%
% Version 1
% Version 2
%
close all
clear
clc

% define constants
R=287.058;       % specific gas constant for air
gam=1.4;         % ratio of specific heat
cv=0.718e3;      % specific heat for air
kappa=24.35e-3;  % conductivity of air
muref=1.716e-5;  % reference viscosity for air
Tref=273.15;     % reference tempereature for air

% set up problem
L=1.0;
numele=1;  % number of elements
numstp=1;  % number of time steps
dt=1.0;     % time step size

inmax=5;     % maximum number of steps in Newton
epsnew=1e-5; % required convergence in Newton
alpnew=1.0;  % under-relaxation in Newton

Q=0.0;         % volumetric heat load
Qfac=200;      % exponential decay factor

rhoini=1.29;   % initial density
velini=0.0;    % initial velocity
tmpini=273.0;  % initial temperature

ctrho=100.0;     % time continuity penalty factor for density 
ctvel=100.0;     % time continuity penalty factor for velocity
cttmp=100.0;     % time continuity penalty factor for temperature

gamma = 10.0;  % Nitsche penalty parameter
beta = 1.0;     % for symmetric / unsymmetric Nitsche

% generate fe mesh
elen=L/numele;
xp=0:elen:L;

% setup solution vector
numnod=numele+1;
numdof=numnod*6;

solvec=zeros(numdof,1); % dofs 1-3 at t1 dofs 4-6 at t2
                        
% prescribed vel dofs at ends of domain
pdof=[];

% define free dofs
fdof=1:numdof;         
fdof(pdof)=[];

%initialize solution
solvec(1:6:end)=rhoini;
solvec(4:6:end)=rhoini;
solvec(2:6:end)=velini;
solvec(5:6:end)=velini;
solvec(3:6:end)=tmpini;
solvec(6:6:end)=tmpini;

% compute nodal heat load
qnodal=zeros(numnod,1);
for in=1:numnod
    x=(in-1)*elen;
    %qnodal(in)=Q*exp(-Qfac*(2*x/L-1)^2);
    qnodal(in)=Q;
end

% limits for plots
rmin=0.8*rhoini;
rmax=1.05*rhoini;

vmin=-2.0;
vmax= 2.0;

tmin=0.9*tmpini;
tmax=2*tmpini;

pmin=1.0*R*rhoini*tmpini;
pmax=1.1*R*rhoini*tmpini;

machmin=0;
machmax=0.01;

reymin=0;
reymax=1e5;

% time loop
for it=1:numstp
    
    % compute time instances
    t1=(it-1)*dt;
    t2=t1+dt;
    t1prev=t1-dt;
    t2prev=t1;
    
    % store solution of previous time slab
    prevec=solvec;
    
    % Newton loop
    for in=1:inmax
        
        res=zeros(numdof,1);
        jac=zeros(numdof,numdof);
        
        % build residual and jacobian
        for ie=1:numele
            
            % get elemental dofs and solution
            edof=(ie-1)*6+[1 7 4 10  2 8 5 11  3 9 6 12];
            esol=solvec(edof);
            eprev=prevec(edof);
                       
            % compute elemetal residual and jacobian from bulk integrals
            [eres,ejac]=ondim_func(esol,0,elen,t1,t2,qnodal(ie),qnodal(ie+1) );
            
            % compute elemetal residual and jacobian from bulk integrals
            [etres,etjac]=ondim_tc_func(esol,eprev, ...
                0,elen,t1,t2,t1prev,t2prev,    ...
                ctrho, ctvel, cttmp );
            
            % add time continuity contribution to residual / jacobian
            eres = eres + etres;
            ejac = ejac + etjac;            
            
            % compute elemental residual and jacobian contribution from
            % Nitsche
            if ( ie == 1 )
                [enres,enjac] = nitsche_func(esol,0,elen,t1,t2,-1,beta,gamma);
                eres = eres + enres;
                ejac = ejac + enjac;
                                                
                % debug
                %[enres,enjac,enjacAna] = nitsche_func(esol,0,elen,t1,t2,-1,beta,gamma);
                %test = enjac - enjacAna
            
            elseif ( ie == numele )
                [enres,enjac] = nitsche_func(esol,0,elen,t1,t2,1,beta,gamma);
                eres = eres + enres;
                ejac = ejac + enjac;
                
                % debug
                %[enres,enjac,enjacAna] = nitsche_func(esol,0,elen,t1,t2,1,beta,gamma);
                %test = enjac - enjacAna
            end
        
            % assemble into global residual and stiffness matrix
            res(edof)=res(edof)+eres;
            jac(edof,edof)=jac(edof,edof)+ejac;
        end
        
        % save reference residual
        if in==1
            refRes=norm(res(fdof));
        end
        
        fprintf('time-step %d new-iter %d  res-norm %e  rel-res-norm %e sol-norm %e cond %e\n', ...
            it,in,norm(res(fdof)),norm(res(fdof))/refRes,norm(solvec),cond(jac(fdof,fdof)));
        
        % check for convergence
        if norm(res(fdof)) < epsnew * refRes || norm(res(fdof)) < 1e-9
            break;
        end
        
        % solve system and update solution
        solvec(fdof)=solvec(fdof)-alpnew*jac(fdof,fdof)\res(fdof);
    end
    
    if in==inmax
        display(['Newton did not converge - relative residual is ' num2str(norm(res(fdof))/refRes)]);
    end
    
    % compute nodal solutions at center for time step
    nrho=zeros(numnod,1);
    nvel=zeros(numnod,1);
    ntmp=zeros(numnod,1);

    drho=0;
    dvel=0;
    dtmp=0;
    
    for in=1:numnod
        idof=(in-1)*6;
        nrho(in)=(solvec(idof+1)+solvec(idof+4))/2;
        nvel(in)=(solvec(idof+2)+solvec(idof+5))/2;
        ntmp(in)=(solvec(idof+3)+solvec(idof+6))/2;
        
        drho=drho+(solvec(idof+1)-solvec(idof+4))^2;
        dvel=dvel+(solvec(idof+2)-solvec(idof+5))^2;
        dtmp=dtmp+(solvec(idof+3)-solvec(idof+6))^2;
    end
    
    fprintf('\ntime-step %d time jump in percent: density %e velocity %e temperature %e\n\n', ...
        it,...
        sqrt(drho)/norm(nrho)*100, ...
        sqrt(dvel)/norm(nvel)*100, ...
        sqrt(dtmp)/norm(ntmp)*100);
 
    
    % compute local Mach and Reynolds numbers
    mach=zeros(numnod,1);
    reyn=zeros(numnod,1);
    pres=zeros(numnod,1);
    for in=1:numnod
        %nmu=mu_func(ntmp(in));     % local viscosity
        nmu = muref;
        ncs=sqrt(gam*R*ntmp(in));  % local speed of sound
        
        reyn(in)=abs(nvel(in))*elen*nrho(in)/nmu; % Reynolds number defined wrt element length
        mach(in)=abs(nvel(in))/ncs;               % Mach number
        pres(in)=R*nrho(in)*ntmp(in);
    end

    h=figure(1);
    
    subplot(3,3,1);
    plot(xp,nrho,'r-');
    axis([0 L rmin rmax]);
    title('density')
    
    subplot(3,3,2);
    plot(xp,nvel,'b-');
    axis([0 L vmin vmax]);
    title('velocity')
    
    subplot(3,3,3);
    plot(xp,ntmp,'k-');
    axis([0 L tmin tmax]);
    title('temperature')

    subplot(3,3,4);
    plot(xp,mach,'r-');
    axis([0 L machmin machmax]);
    title('local Mach number')

    subplot(3,3,5);
    plot(xp,reyn,'r-');
    axis([0 L reymin reymax]);
    title('local Reynolds number')
    
    subplot(3,3,6);
    plot(xp,pres,'g-');
    axis([0 L pmin pmax]);
    title('pressure')
    
    subplot(3,3,7);
    plot(xp,qnodal,'g-');
    title('nodal heat flux')
    
    getframe(h);
    
end


