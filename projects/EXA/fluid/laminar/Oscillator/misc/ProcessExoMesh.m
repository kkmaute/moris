function [MassMat,fecoord,fetopo,feprop,febc,numfenod,numfele]=ProcessExoMesh(filename,plotflag)
%
% matlab script to read mesh and compute basic fe-data:
% such as mass matrix, volume, area .....
%
% input:   filename - name of exo file (only 4-node elements possible)
%          plotflag - flag for plotting mesh:
%                     < 1  - no plot
%                       1  - just mesh
%                       2  - mehs with node and element numbers
%
% output:  MassMat  - mass matrix in sparse format
% 

% set path (need to be run in Transfem root directory)
addpath(genpath(pwd));

% hardwired parameters

ndnc=1;    % number of dofs per node

% read mesh

[fecoord,fetopo,feprop,febc,numfenod,numfele]=LoadExoMesh(filename);

% generate dofs for 4-node 

ndex=4*ndnc;                            % number of nodal dofs per element
numdof=ndnc*numfenod;                   % total number of nodal dofs in mesh (will be updated later)
ndof=reshape(1:numdof,ndnc,numfenod)';  % list of nodal dof ids

ex=zeros(numfele,4);
ey=zeros(numfele,4);

edof=zeros(numfele,ndex+1);
edof(:,1)=1:numfele;

% set nodal degree of freedom

for i=1:numfele
     ex(i,:)=fecoord(fetopo(i,2:5),1);
     ey(i,:)=fecoord(fetopo(i,2:5),2);
     edof(i,2:ndex+1)=reshape(ndof(fetopo(i,2:5),:)',1,ndex); 
end

% plot mesh
if plotflag > 0
    eldraw2(ex,ey,[1 1 0]);
    if plotflag > 1
        for ii=1:numfenod
            h=text(fecoord(ii,1),fecoord(ii,2),int2str(ii));
            set(h,'fontsize',10,'color','r');
        end
        for ii=1:numfele
            h=text(mean(ex(ii,:)),mean(ey(ii,:)),int2str(ii));
            set(h,'fontsize',10,'color','b');
        end
    end
end

% build fe-data

ndele=4*ndnc;
ndele2=ndele^2;
nzn=numfele*ndele;

MassMatRow=zeros(nzn,1);
MassMatCol=zeros(nzn,1);
MassMatVec=zeros(nzn,1);

nnz=0;

for iele=1:numfele
    
    [me]=flw2i4m(ex(iele,:),ey(iele,:),[1 2],1,1);
    
    ind1       = (nnz+1):(nnz+ndele2);
    ind2       = edof(iele,2:ndele+1);
    ind3       = repmat(ind2,ndele,1);
    
    nnz        = nnz+ndele2;
    
    MassMatCol(ind1) = reshape(ind3,  1,ndele2);
    MassMatRow(ind1) = reshape(ind3', 1,ndele2);
    MassMatVec(ind1) = reshape(me,    1,ndele2);
    
end

MassMat = sparse(MassMatRow,MassMatCol,MassMatVec,numdof,numdof,nnz);
