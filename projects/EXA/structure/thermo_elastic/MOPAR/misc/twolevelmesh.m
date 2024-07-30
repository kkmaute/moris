function [fecoord_fine,topo_fine,fecoord_crse,topo_crse,Mmat,bmat]= ...
    twolevelmesh(xlen,ylen,nx_fine,ny_fine,nx_crse,ny_crse)
%==========================================================================
% 
% 
%==========================================================================

% compute mesh parameters
dx_fine=xlen/nx_fine;   % fine mesh edge length (x-direction)
dy_fine=ylen/ny_fine;   % fine mesh edge length (y-direction)

dx_crse=xlen/nx_crse;   % coarse mesh edge length (x-direction)
dy_crse=ylen/ny_crse;   % coarse mesh edge length (y-direction)

% create nodal coordinates
numnod_fine=(nx_fine+1)*(ny_fine+1); % number of nodes in fine mesh
numnod_crse=(nx_crse+1)*(ny_crse+1); % number of nodes in coarse mesh

fecoord_fine=zeros(numnod_fine,2);   % nodal coordinates in fine mesh
fecoord_crse=zeros(numnod_crse,2);   % nodal coordinates in coarse mesh

nz=0;
for iy=1:ny_fine+1
    for ix=1:nx_fine+1
        nz=nz+1;
        fecoord_fine(nz,:)=[(ix-1)*dx_fine (iy-1)*dy_fine];
    end
end

nz=0;
for iy=1:ny_crse+1
    for ix=1:nx_crse+1
        nz=nz+1;
        fecoord_crse(nz,:)=[(ix-1)*dx_crse (iy-1)*dy_crse];
    end
end

% create elements
numele_fine=nx_fine*ny_fine;  % number of elements in fine mesh
numele_crse=nx_crse*ny_crse;  % number of elements in coarse mesh

topo_fine=zeros(numele_fine,4+1);
topo_crse=zeros(numele_crse,4+1);

topo_fine(:,1)=3;
topo_crse(:,1)=3;

nz=0;
for iy=1:ny_fine
    for ix=1:nx_fine
        nz=nz+1;
        baseId=(nx_fine+1)*(iy-1)+ix;
        topo_fine(nz,2:end)=baseId+[ 0    1       nx_fine+2      nx_fine+1];
    end
end

nz=0;
for iy=1:ny_crse
    for ix=1:nx_crse
        nz=nz+1;
        baseId=(nx_crse+1)*(iy-1)+ix;
        topo_crse(nz,2:end)=baseId+[ 0    1       nx_crse+2      nx_crse+1];
    end
end

% create mapping operators assuming bi-linear interpolations
mm_row=zeros(numele_fine*16*4,1);
mm_col=zeros(numele_fine*16*4,1);
mm_val=zeros(numele_fine*16*4,1);

bb_row=zeros(numele_fine*16*4,1);
bb_col=zeros(numele_fine*16*4,1);
bb_val=zeros(numele_fine*16*4,1);

nz=0;
pz=0;
% loop over all fine elements
for ie=1:numele_fine
    
    % node ids of element (fine mesh)
    noids_fine = topo_fine(ie,2:end);
    
    % get elemental coordinates (fine mesh)
    ex_fine=fecoord_fine(noids_fine,1);
    ey_fine=fecoord_fine(noids_fine,2);
    
    % locations get gauss points (global and local) and weights in fine mesh
    [xgp_fine,ygp_fine,wgp_fine]=xfem4quad(ex_fine',ey_fine',0,0);
    
    % loop over all gauss points
    for ig=1:length(wgp_fine)
        
        % get shape functions in fine mesh
        [N_fine]=xfem4stdshape(xgp_fine(ig,:),ygp_fine(ig,:));
        
        % determine element id in coarse mesh
        ix_crse=floor(xgp_fine(ig,1)/dx_crse);
        iy_crse=floor(ygp_fine(ig,1)/dy_crse);
        
        ie_crse=nx_crse*iy_crse+ix_crse+1;
        
        % node ids of element (corse mesh)
        noids_crse = topo_crse(ie_crse,2:end);
        
        % get elemental coordinates (coarse mesh)
        ex_crse=fecoord_crse(noids_crse,1);
        ey_crse=fecoord_crse(noids_crse,2);

        % find local coordinates in corse element 
        % (assume rectangular element)
        xgp_crse=(2*xgp_fine(ig,1)-ex_crse(1)-ex_crse(2))/(ex_crse(2)-ex_crse(1));
        ygp_crse=(2*ygp_fine(ig,1)-ey_crse(2)-ey_crse(3))/(ey_crse(3)-ey_crse(2));
        
        if abs(xgp_crse) > 1 || abs(ygp_crse) > 1; 
            error('outslide bounds'); 
        end
        
        % get shape functions in coarse mesh
        [N_crse]=xfem4stdshape([0 xgp_crse],[0 ygp_crse]);
        
        % integrate error function
        mm = N_fine'*N_fine*wgp_fine(ig);
        bb = N_fine'*N_crse*wgp_fine(ig);
        
        % store entries in sparse matrix form
        ind_fine = repmat(noids_fine,4,1);
        mm_row(nz+1:nz+16)=reshape(ind_fine',1,16);
        mm_col(nz+1:nz+16)=reshape(ind_fine ,1,16);
        mm_val(nz+1:nz+16)=reshape(mm,1,16);
        nz=nz+16;
        
        ind_crse = repmat(noids_crse,4,1);
        bb_row(pz+1:pz+16)=reshape(ind_fine',1,16);
        bb_col(pz+1:pz+16)=reshape(ind_crse ,1,16);
        bb_val(pz+1:pz+16)=reshape(bb,1,16);
        pz=pz+16;
    end
end
   
% build sparse matrices
Mmat=sparse(mm_row(1:nz),mm_col(1:nz),mm_val(1:nz),numnod_fine,numnod_fine,nz);
bmat=sparse(bb_row(1:pz),bb_col(1:pz),bb_val(1:pz),numnod_fine,numnod_crse,pz);

end
