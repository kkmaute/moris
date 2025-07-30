% 
%   Copyright (c) 2023 University of Colorado
%   Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%  
%  ------------------------------------------------------------------------------------
%==========================================================================
% script to extract edege from multi-block set exodus file and write 
% it to 2D obj file
%==========================================================================
function ExtractEdges

clc
close all
clear all

%==========================================================================
%
% User defined parameters
%
%==========================================================================

% base name of exodus output files
outputbase='FSC.exo.e-s.0501';

% number of processors
nproc=16;

% phase to be extracted
phase=1;

% tolerance for merging nodes
deps=1e-6;

% target length for resampleing segments 
newseg=6e-1;

% smoothing
smonum=25;
smostep=[0.33,-0.34];

% flag for computing inscribed circle diameter
iscflg=1;

% flag for ploting structure
pltflg=0;

% flag for ploting inscribed circles
pltcirc=0;

% file of output obj file
objfilename='extracted.obj';

%==========================================================================
%
% Endo fo User defined parameters
%
%==========================================================================

% read mesh and 
[crds,quads,tris]=readmesh(outputbase,nproc,phase);

% remove coinciding nodes
[crds,quads,tris]=cleanmesh(crds,quads,tris,deps);

% extract number of nodes, quads, and tris
numnodes=size(crds,1);
numquads=size(quads,1);
numtris=size(tris,1);

% plot mesh
if pltflg>0
    figure(1)
    for i=1:numquads
        nids=[quads(i,:) quads(i,1)];
        xp=crds(nids,1);
        yp=crds(nids,2);
        plot(xp,yp,'k-'); hold on;
    end
    for i=1:numtris
        nids=[tris(i,:) tris(i,1)];
        xp=crds(nids,1);
        yp=crds(nids,2);
        plot(xp,yp,'k-'); hold on;
    end
    axis equal
end

% extract edges
maxedg=4*numquads+3*numtris;
edgcon=zeros(maxedg,3);
edgconorg=zeros(maxedg,2);
ic=0;
for i=1:numquads
    nids=[quads(i,:) quads(i,1)];
    for j=1:4
        ni=quads(i,j);
        nj=nids(j+1);
        ic=ic+1;
        edgcon(ic,:)=[sort([ni nj],'ascend') i];
        edgconorg(ic,:)=[ni nj];
    end
end
for i=1:numtris
    nids=[tris(i,:) tris(i,1)];
    for j=1:3
        ni=tris(i,j);
        nj=nids(j+1);
        ic=ic+1;
        edgcon(ic,:)=[sort([ni nj],'ascend') i+numquads];
        edgconorg(ic,:)=[ni nj];
    end
end

% identify edges connected only to a single element
iec=0;
exedge=[];
unodes=zeros(numnodes,1);
for i=1:size(edgcon,1)
    ni=edgcon(i,1);
    nj=edgcon(i,2);
    single=1;
    cedg=find(edgcon(:,1)==ni);
    for j=1:length(cedg)
        ig=cedg(j);
        njj=edgcon(ig,2);
        if njj ~= nj; continue; end
        if edgcon(i,3) ~= edgcon(ig,3)
            single=0;
        end
    end
    if single==1
        iec=iec+1;
        ni=edgconorg(i,1);
        nj=edgconorg(i,2);
        unodes([ni nj])=1;
        exedge=[exedge; ni nj];
    end
end

fprintf('\n ... number of extracted edges: %d\n',iec);

% build edge mesh
inc=0;
excrd=[];
for i=1:numnodes
    if unodes(i)>0
        inc=inc+1;
        unodes(i)=inc;
        excrd=[excrd;crds(i,:)];
    end
end

for i=1:iec
    exedge(i,:)=unodes(exedge(i,:));
end

if pltflg>0
    figure(2)
    for i=1:iec
        ni=exedge(i,1);
        nj=exedge(i,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);
        plot(xp,yp,'k-'); hold on;
    end
    axis equal
end

% resample segments
if newseg>0
     [excrd,exedge]=removesegements(excrd,exedge,newseg,pltflg);
end

% smooth shape
if smonum>0
     excrdorg=excrd;
     [excrd,dexcrd,freezenodes]=smoothing(excrd,exedge,smostep,smonum,pltflg);
     
     numnod=size(excrd,1);
     fderror=zeros(numnod,2);
     deps=1e-5;
     maxerror=0;
     for i=1:numnod
         if freezenodes(i);continue,end
         for j=1:2
             nper=i;
             cper=j;

             excrdp=excrdorg;
             excrdp(nper,cper)=excrdp(nper,cper)+deps;
             [excrdp]=smoothing(excrdp,exedge,smostep,smonum,0);

             excrdm=excrdorg;
             excrdm(nper,cper)=excrdm(nper,cper)-deps;
             [excrdm]=smoothing(excrdm,exedge,smostep,smonum,0);
             fd=(excrdp-excrdm)/2/deps;
             fderror(i,j)=norm(dexcrd(:,nper)-fd(:,cper))/norm(fd(:,cper))*100.0;
             if maxerror<fderror(i,j)
                 maxerror=fderror(i,j);
                 nodmax=i;
                 crdmax=j;
             end
         end
     end
     fprintf('\n ... max smoothing sensitivity error in percent: %e (node %d dir %d)\n',...
         maxerror,nodmax,crdmax);
end

% compute inscribed circle diameter
if iscflg>0
    [iscdiameter]=computeInscribedCircleDiameter(excrd,exedge,pltflg,pltcirc);
end

% write extracted edges to file
fid = fopen(objfilename,'w');

for i = 1:size(excrd,1)
    fprintf(fid, 'v %f %f\n', excrd(i,1), excrd(i,2));
end
for i = 1:size(exedge,1)
    fprintf(fid, 'f %d %d\n', exedge(i,1), exedge(i,2));
end
fclose(fid);

fprintf('\n ... created obj file: %s\n',objfilename);

end

%==========================================================================

function [crds,quads,tris]=readmesh(outputbase,nproc,phase)

% initialize arrays
quads=[];
tris=[];
fecs=[];
ndmap=[];

% loop over all processors
for ip=1:nproc

    if nproc>1
        if nproc >= 100
            outputmesh=sprintf('%s.%d.%03d',outputbase,nproc,ip-1);
        else
            if nproc >= 10
                outputmesh=sprintf('%s.%d.%02d',outputbase,nproc,ip-1);
            else
                outputmesh=sprintf('%s.%d.%d',outputbase,nproc,ip-1);
            end
        end
    else
        outputmesh=sprintf('%s.exo',outputbase);
    end

    fprintf(' ... processing file: %s\n',outputmesh);

    % extract mesh information from exodus file
    MeshData = ExtractExoData(0,outputmesh);

    % nodal coordiantes
    fecoord=zeros(length(MeshData.Coord{1,1}.Coords),length(MeshData.Coord));
    for i=1:1:length(MeshData.Coord)
        fecoord(:,i) = MeshData.Coord{1,i}.Coords;
    end

    fecs=[fecs; fecoord];

    % node map
    ndm=MeshData.NodeMap+1;
    ndmap=[ndmap; ndm];

    % number of element blocks
    numblk=size(MeshData.Etopo,2);

    for i=1:numblk
        bsname=MeshData.Etopo{1,i}.Name;
        pid=str2num(bsname(14));
        if pid==phase
            npe=size(MeshData.Etopo{1,i}.Nodes,1);
            if npe>3
                quads=[quads; ndm(MeshData.Etopo{1,i}.Nodes(1:4,:))'];
            else
                tris=[tris; ndm(MeshData.Etopo{1,i}.Nodes(1:3,:))'];
            end
        end
    end
end

% extract number of nodes
numnodes=length(ndmap);

% get maximum node id
mxndid=max(ndmap);
crds=zeros(mxndid,2);

% build nodal coordinates with global node ids
for i=1:numnodes
    nid=ndmap(i);
    crds(nid,:)=fecs(i,:);
end

end

%==========================================================================

function [crds,quads,tris]=cleanmesh(crds,quads,tris,deps)

% number of bins in direction with maximum dimension
nx=200;

% extract number of nodes, quads, and tris
numnodes=size(crds,1);
numquads=size(quads,1);
numtris=size(tris,1);

% sort nodes into background boxes
xmin=min(crds(:,1));
ymin=min(crds(:,2));
xdist=max(crds(:,1))-xmin;
ydist=max(crds(:,2))-ymin;
dx=max(xdist,ydist)/nx;
deps=deps*max(xdist,ydist);
nx=ceil(xdist/dx);
ny=ceil(ydist/dx);
dx=xdist/nx;
dy=ydist/ny;

% count number of points in each box
npib=zeros(nx,ny);
for i=1:numnodes
    ix=min(floor((crds(i,1)-xmin)/dx)+1,nx);
    iy=min(floor((crds(i,2)-ymin)/dy)+1,ny);
    npib(ix,iy)=npib(ix,iy)+1;

    %add to sourrounding boxes
    for m=-1:1
        for n=-1:1
            if m==0 && n==0; continue;end
            ixx=ix+m;
            iyy=iy+n;
            if ixx>=1 && ixx<=nx && iyy>=1 && iyy<=ny
                npib(ixx,iyy)=npib(ixx,iyy)+1;
            end
        end
    end
end

nibx=cell(nx,ny);
for i=1:nx
    for j=1:ny
        nibx{i,j}=zeros(npib(i,j),1);
    end
end

npib=zeros(nx,ny);
nbpn=zeros(numnodes,1);
for i=1:numnodes
    ix=min(floor((crds(i,1)-xmin)/dx)+1,nx);
    iy=min(floor((crds(i,2)-ymin)/dy)+1,ny);
    npib(ix,iy)=npib(ix,iy)+1;
    nibx{ix,iy}(npib(ix,iy))=i;
    nbpn(i)=nbpn(i)+1;

    %add to sourrounding boxes
    for m=-1:1
        for n=-1:1
            if m==0 && n==0; continue;end
            ixx=ix+m;
            iyy=iy+n;
            if ixx>=1 && ixx<=nx && iyy>=1 && iyy<=ny
                npib(ixx,iyy)=npib(ixx,iyy)+1;
                nibx{ixx,iyy}(npib(ixx,iyy))=i;
                nbpn(i)=nbpn(i)+1;
            end
        end
    end
end

bfen=cell(numnodes,1);
for i=1:numnodes
    bfen{i}=zeros(nbpn(i),2);
end

nbpn=zeros(numnodes,1);
for i=1:nx
    for j=1:ny
        for k=1:length(nibx{i,j})
            nid=nibx{i,j}(k);
            nbpn(nid)=nbpn(nid)+1;
            bfen{nid}(nbpn(nid),:)=[i j];
        end
    end
end

% find and remove coinciding nodes
ndmap=zeros(numnodes,1);
ic=0;
for i=1:numnodes
    if ndmap(i) > 0; continue;end
    ic=ic+1;
    ndmap(i)=ic;

    for ib=1:nbpn(i)
        ix=bfen{i}(ib,1);
        iy=bfen{i}(ib,2);
        nodlist=nibx{ix,iy};
        for j=1:length(nodlist)
            nj=nodlist(j);
            if ndmap(nj) > 0; continue;end
            dist=norm(crds(i,:)-crds(nj,:));
            if dist<deps
                ndmap(nj)=ic;
            end
        end
    end
end

fprintf('\n ... number of redundant nodes: %d\n',numnodes-ic);

newcrds=zeros(ic,2);
for i=1:numnodes
    nnid=ndmap(i);
    newcrds(nnid,:)=crds(i,:);
end
crds=newcrds;
for i=1:numquads
    quads(i,:)=ndmap(quads(i,:));
end
for i=1:numtris
    tris(i,:)=ndmap(tris(i,:));
end
end

%==========================================================================

function Data = ExtractExoData(printflag,datafile,savefile)

% open the specified netCDF file or reading
ncid = netcdf.open(datafile,'NC_NOWRITE');

% Get the number of dimensions and variables saved in the file
[ndims,nvars] = netcdf.inq(ncid);

% Scan through and print all dimension names and their values
Dim{1}.Name = '';
Dim{1}.Value = '';
for id = 0:ndims-1
    [Dim{id+1}.Name, Dim{id+1}.Value] = netcdf.inqDim(ncid,id);
    if printflag == 1
        fprintf(['dimension %d ' Dim{id+1}.Name ' = %d \n'],id,Dim{id+1}.Value)
    end
end

% get number of element blocks and elements in each block
num_el_blk=Dim{6}.Value;
num_el_in_blk=zeros(num_el_blk,1);
for id = 1:ndims
    if length(Dim{id}.Name) == 14
        if strcmp(Dim{id}.Name(1:13),'num_el_in_blk')
            blkid=str2num(Dim{id}.Name(14));
            num_el_in_blk(blkid)=Dim{id}.Value;
        end
    end
end

% initialize variables of interest
ci  = 0;
toi = 0;
nsi = 0;
ssi = 0;
nvi = 0;
evi = 0;
bvi = 0;
nsID_counter = 0;
ssID_counter = 0;
Coordnames = '';
BSnames = '';
NSnames = '';
SSnames = '';
NVnames = '';
EVnames = '';
GVnames = '';
Time    = [];
Coord   = [];
Etopo   = [];
NodeSet = [];
nsID    = [];
SideSet = [];
ssID    = [];
GVar    = [];
NVar    = [];
EVar    = [];
elem_map= [];
node_map= [];
% Scan through and all variables printing their names the names and
% values of any associated attributes, and save the desired variables.
for iv = 0:nvars-1
    varflag = 1;
    varstring = '';
    Att = [];
    % Get the variable name and the number of associated attributes
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,iv);
    % Put spaces on the end of the string to prevent strcmp errors
    text = [varname '                '];

    if printflag == 1
        % Print the variable number and name
        fprintf(['variable %d is ' text '\n'],iv)
    end
    % Print any attribute names and values associated with the current variable
    for ia = 0:natts-1
        Att{ia+1}.Number = ia;
        Att{ia+1}.Name   = netcdf.inqAttName(ncid,iv,ia);
        Att{ia+1}.Value  = netcdf.getAtt(ncid,iv,netcdf.inqAttName(ncid,iv,ia));
        if printflag == 1
            fprintf(['\tattribute %d ' netcdf.inqAttName(ncid,iv,ia) ' = ' netcdf.getAtt(ncid,iv,netcdf.inqAttName(ncid,iv,ia)) '\n'],ia)
        end
    end

    % Save desired variables
    if strcmp(text(1:10),'time_whole')
        Time = netcdf.getVar(ncid,iv);
        varstring = 'Time';
    elseif strcmp(text(1:5),'coord')
        coordinates = netcdf.getVar(ncid,iv);
        dcidim = size(coordinates,2);
        for idc = 1:dcidim
            ci = ci+1;
            Coord{ci}.Coords = coordinates(:,idc);
        end
        varstring = sprintf('Coord{}.Coords');
    elseif strcmp(text(1:10),'coor_names')
        Coordnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('Coord{}.Name');
    elseif strcmp(text(1:7),'connect')
        toi = toi+1;
        for ia = 0:natts-1
            attname = netcdf.inqAttName(ncid,iv,ia);
            text1 = [attname '          '];
            if strcmp(text1(1:9),'elem_type')
                Etopo{toi}.EleType = netcdf.getAtt(ncid,iv,netcdf.inqAttName(ncid,iv,ia));
            end
        end
        Etopo{toi}.Nodes = netcdf.getVar(ncid,iv);
        varstring = sprintf('Etopo{%d}.Nodes',toi);
    elseif strcmp(text(1:7),'node_ns')
        nsi = nsi+1;
        NodeSet{nsi}.Nodes = netcdf.getVar(ncid,iv);
        varstring = sprintf('NodeSet{%d}.Nodes',nsi);
    elseif strcmp(text(1:12),'dist_fact_ns')
        NodeSet{nsi}.DistFact = netcdf.getVar(ncid,iv);
        varstring = sprintf('NodeSet{%d}.DistFact',nsi);
    elseif strcmp(text(1:8),'eb_names')
        BSnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('Etopo{}.Name');
    elseif strcmp(text(1:8),'ns_names')
        NSnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('NodeSet{}.Name');
    elseif strcmp(text(1:7),'ns_prop')
        nsID_counter = nsID_counter +1;
        nsID{nsID_counter} = netcdf.getVar(ncid,iv);
    elseif strcmp(text(1:7),'elem_ss')
        ssi = ssi+1;
        SideSet{ssi}.Elements = netcdf.getVar(ncid,iv);
        varstring = sprintf('SideSet{%d}.Elements',ssi);
    elseif strcmp(text(1:7),'side_ss')
        SideSet{ssi}.Sides = netcdf.getVar(ncid,iv);
        varstring = sprintf('SideSet{%d}.Sides',ssi);
    elseif strcmp(text(1:7),'ss_prop')
        ssID_counter = ssID_counter +1;
        ssID{ssID_counter} = netcdf.getVar(ncid,iv);
    elseif strcmp(text(1:12),'dist_fact_ss')
        SideSet{ssi}.DistFact = netcdf.getVar(ncid,iv);
        varstring = sprintf('SideSet{%d}.DistFact',ssi);
    elseif strcmp(text(1:8),'ss_names')
        SSnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('SideSet{}.Name');
    elseif strcmp(text(1:12),'vals_glo_var')
        GVarVal = netcdf.getVar(ncid,iv);
        varstring = sprintf('GVar{}.Value');
    elseif strcmp(text(1:12),'name_glo_var')
        GVnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('GVar{}.Name');
    elseif strcmp(text(1:12),'vals_nod_var')
        nvi = nvi+1;
        NVar{nvi}.Val = netcdf.getVar(ncid,iv);
        varstring = sprintf('NVar{%d}.Value',nvi);
    elseif strcmp(text(1:12),'name_nod_var')
        NVnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('NVar{}.Name');
    elseif strcmp(text(1:13),'vals_elem_var')
        kk=strfind(text,'eb');
        varid=str2num(text(14:kk-1));
        blkid=str2num(text(kk+2:end));
        evi = max(evi,varid);
        bvi = max(bvi,blkid);
        EVar{varid,blkid}.Val = netcdf.getVar(ncid,iv);
        varstring = sprintf('EVar{%d}.Value',varid);
    elseif strcmp(text(1:13),'name_elem_var')
        EVnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('EVar{}.Name');
    elseif strcmp(text(1:12),'elem_num_map')
        elem_map=netcdf.getVar(ncid,iv);
    elseif strcmp(text(1:12),'node_num_map')
        node_map=netcdf.getVar(ncid,iv);
    else
        Var{iv+1}.Number = iv;
        Var{iv+1}.Name   = varname;
        Var{iv+1}.Value  = netcdf.getVar(ncid,iv);
        varflag = 0;
    end
    if varflag
        Var{iv+1}.Number = iv;
        Var{iv+1}.Name   = varname;
        Var{iv+1}.Value  = varstring;
    end
    Var{iv+1}.Att = Att;
end
if ~isempty(Coordnames)
    for ic = 1:2 %ci
        Coord{ic}.Name = Coordnames(:,ic)';
    end
end
if ~isempty(BSnames)
    ic=0;
    for in = 1:num_el_blk
        if num_el_in_blk(in)>0
            ic=ic+1;
            Etopo{ic}.Name = BSnames(:,in)';
        end
    end
end
if ~isempty(NSnames)
    for in = 1:nsi
        NodeSet{in}.Name = NSnames(:,in)';
    end
end
if ~isempty(SSnames)
    for is = 1:ssi
        SideSet{is}.Name = SSnames(:,is)';
    end
end
if ~isempty(GVnames)
    for ig = 1:size(GVarVal,1)
        GVar{ig}.Value = GVarVal(ig,:);
        GVar{ig}.Name = GVnames(:,ig)';
    end
end
if ~isempty(NVnames)
    for in = 1:nvi
        NVar{in}.Name = NVnames(:,in)';
    end
end
if ~isempty(EVnames)
    for in = 1:evi
        for ib = 1:bvi
            EVar{in,ib}.Name = EVnames(:,in)';
        end
    end
end

% close the specified netCDF file
netcdf.close(ncid)

% Create a cell type data structure to be saved
Data.Dim       = Dim;
Data.Var       = Var;
Data.Time      = Time;
Data.Coord     = Coord;
Data.Etopo     = Etopo;
Data.NodeSet   = NodeSet;
Data.NodeSetID = nsID;
Data.SideSet   = SideSet;
Data.SideSetID = ssID;
Data.GVar      = GVar;
Data.NVar      = NVar;
Data.EVar      = EVar;
Data.EleMap    = elem_map;
Data.NodeMap   = node_map;

% if specified save the data in a matlab file
if nargin == 3
    save(savefile,'Data')
end
end

%==========================================================================

function [iscdiameter]=computeInscribedCircleDiameter(points,connect,pltflg,pltcirc)

numnod=size(points,1);
numseg=size(connect,1);

if pltflg>0
    figure(5)
    for i=1:numseg
        ni=connect(i,1);
        nj=connect(i,2);
        xp=points([ni nj],1);
        yp=points([ni nj],2);
        plot(xp,yp,'k-'); hold on;
    end
    axis equal
end

% get bounding box
xmin=min(points(:,1));
xmax=max(points(:,1));
ymin=min(points(:,2));
ymax=max(points(:,2));

maxradius=sqrt((xmax-xmin)^2+(ymin-ymax)^2);

% compute nodal connectivity
nodconnect=zeros(numnod,2);
for i=1:numseg
    nodeA = connect(i, 1);
    nodeB = connect(i, 2);
    if nodconnect(nodeA,1)==0; k=1; else; k=2; end
    nodconnect(nodeA,k) = i;
    if nodconnect(nodeB,1)==0; k=1; else; k=2; end
    nodconnect(nodeB,k) = i;
end

% compute normal on segments
segnormal=zeros(numseg,2);
seglength=zeros(numseg,1);
for i=1:numseg
    nodeA=connect(i,1);
    nodeB=connect(i,2);

    vABx=points(nodeB,1)-points(nodeA,1);
    vABy=points(nodeB,2)-points(nodeA,2);

    % Calculate the length of the segment
    seglength(i) = sqrt(vABx^2 + vABy^2);

    % Compute the normal vector
    segnormal(i,:) = [-vABy, vABx] / seglength(i);
end

% compute nodal normal
nodnormal=zeros(numnod,2);
for i=1:numnod
    segA=nodconnect(i,1);
    segB=nodconnect(i,2);
    nodnormal(i,1)=(segnormal(segA,1)*seglength(segA)+segnormal(segB,1)*seglength(segB))/ ...
        (seglength(segA)+seglength(segB));
    nodnormal(i,2)=(segnormal(segA,2)*seglength(segA)+segnormal(segB,2)*seglength(segB))/ ...
        (seglength(segA)+seglength(segB));

    nodnormal(i,:)=1/norm(nodnormal(i,:))*nodnormal(i,:);
end

% compute inscribed circel
fprintf('\n ... start computing inscribed circles\n')
tic
mrad=zeros(numnod,1);
for i=1:numnod
    minradius=maxradius;
    for j=1:numnod

        % skip node
        if i==j; continue; end

        % compute delta x and y
        deltax=points(i,1)-points(j,1);
        deltay=points(i,2)-points(j,2);

        radius=-(deltax^2+deltay^2)/ ...
            (2*(nodnormal(i,1)*deltax+nodnormal(i,2)*deltay));

        if radius<minradius && radius>0
            minradius = radius;
        end
    end

    if minradius>0.9*maxradius || pltcirc
        plotcircle(points(i,:),nodnormal(i,:),minradius);
    end
    
    mrad(i)=minradius;
end
toc

iscdiameter=2*mrad;

if pltflg
    figure(6)
    scatter(points(:,1), points(:,2), 50, iscdiameter, 'filled');
    colorbar;
    xlabel('X');
    ylabel('Y');
    title('Points Colored by diameter of inscribed circle');
    axis equal
end

end

%=========================================================================
function plotcircle(point,normal,radius)

figure(5)

centerX=point(1)+radius*normal(1);
centerY=point(2)+radius*normal(2);

theta = linspace(0, 2*pi, 100); % Parameter for circle
circleX = radius * cos(theta) + centerX; % X coordinates of the circle
circleY = radius * sin(theta) + centerY; % Y coordinates of the circle
plot(circleX, circleY, 'r-'); hold on;% Plot the inscribed circle
axis equal

end

%==========================================================================

function [excrd,exedge]=removesegements(excrd,exedge,newseg,pltflg)

numnod=size(excrd,1);
numseg=size(exedge,1);

% compute nodal connectivity
nodconnect=zeros(numnod,2);
for i=1:numseg
    nodeA = exedge(i, 1);
    nodeB = exedge(i, 2);
    if nodconnect(nodeA,1)==0; k=1; else; k=2; end
    nodconnect(nodeA,k) = i;
    if nodconnect(nodeB,1)==0; k=1; else; k=2; end
    nodconnect(nodeB,k) = i;
end

% identify closed polygon
segcounter=0;
segid=zeros(numseg,1);

nonassigned=find(segid==0);

% loop over all closed polygons
maxpoly=100;
segcon=cell(maxpoly,1);

while length(nonassigned)>0

    segcounter=segcounter+1;
    polycounter=1;

    segcon(segcounter)=[];
    
    % food fill
    currentseg=nonassigned(1);
    
    while currentseg>0
        segid(currentseg)=segcounter;
        
        segcon{segcounter}(polycounter)=currentseg;
        polycounter=polycounter+1;

        node2=exedge(currentseg,2);
        segs=nodconnect(node2,:);

        currentseg=0;
        if segid(segs(1))==0; currentseg=segs(1);end
        if segid(segs(2))==0; currentseg=segs(2);end
    end

    nonassigned=find(segid==0);
end

if pltflg>0
    figure(3)
    for i=1:numseg
        ni=exedge(i,1);
        nj=exedge(i,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);

        surface([xp xp], [yp yp],zeros(2,2),segid(i)*ones(2,2), ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2); hold on
    end
    axis equal
    colorbar;
end

% resample segment
edgcounter=0;
nodcounter=0;

newcrd=[];
newedge=[];
newsegid=[];

for ipoly=1:segcounter

    numseg=length(segcon{ipoly});

    % compute total segment length
    tlen=0;
    for iseg=1:numseg
        sid=segcon{ipoly}(iseg);
        ni=exedge(sid,1);
        nj=exedge(sid,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);
        tlen=tlen+sqrt((xp(2)-xp(1))^2+(yp(2)-yp(1))^2);
    end

    % skip polygon if too short
    if floor(tlen/newseg)<2; continue; end

    % compute adjusted new segment length
    adjseg=tlen/floor(tlen/newseg);

    % first point on polygon
    cursid=segcon{ipoly}(1);
    curnid=exedge(cursid,1);

    nodcounter=nodcounter+1;
    newcrd(nodcounter,:)=excrd(curnid,:);

    startnode=nodcounter;

    % serach for new node
    iseg=1;
    sold=0;
    while iseg<=numseg

        found=0;

        % find intersction point on current segment
        sid=segcon{ipoly}(iseg);
        ni=exedge(sid,1);
        nj=exedge(sid,2);
        dxab=excrd(nj,:)-excrd(ni,:);
        dxrt=excrd(ni,:)-newcrd(nodcounter,:);

        % solve quadratic equation
        c=dxrt*dxrt'-adjseg^2;
        b=2*dxab*dxrt';
        a=dxab*dxab';
        s=(-b+sqrt(b^2-4*a*c))/2/a;
        if s>=0 && s<=1 && s>sold
            sold=s;
            found=1;
        else
            s=(-b-sqrt(b^2-4*a*c))/2/a;
            if s>=0 && s<=1 && s>sold
                found=1;
                sold=s;
            end
        end

        % register new point and create new segment
        if found==1
            nodcounter=nodcounter+1;
            newcrd(nodcounter,:)=excrd(ni,:)+s*dxab;

            edgcounter=edgcounter+1;
            newedge(edgcounter,:)=[nodcounter-1 nodcounter];
            newsegid(edgcounter)=ipoly;
        else
            iseg=iseg+1;
            sold=0;
        end
    end

    % check if polygon is close
    if norm(newcrd(startnode,:)-newcrd(nodcounter,:)) < 1e-6*adjseg
        nodcounter=nodcounter-1;
        newedge(edgcounter,:)=[nodcounter startnode];
    else
        edgcounter=edgcounter+1;
        newedge(edgcounter,:)=[nodcounter startnode];
        newsegid(edgcounter)=ipoly;
    end
end

exedge=newedge(1:edgcounter,:);
segid=newsegid(1:edgcounter);
excrd=newcrd(1:nodcounter,:);

fprintf('\n ... number of edges after resampling: %d\n',size(exedge,1));

if pltflg>0
    figure(4)
    for i=1:size(exedge,1)
        ni=exedge(i,1);
        nj=exedge(i,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);

        surface([xp xp], [yp yp],zeros(2,2),segid(i)*ones(2,2), ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2); hold on
    end
    axis equal
    colorbar;
end
end

%==========================================================================

function [excrd,dexcrd,freezenodes]=smoothing(excrd,exedge,smostep,smonum,pltflg)

numnod=size(excrd,1);
numseg=size(exedge,1);

% get bounding box
xmin=min(excrd(:,1));
xmax=max(excrd(:,1));
ymin=min(excrd(:,2));
ymax=max(excrd(:,2));

tolerance=1e-5*max(xmax-xmin,ymax-ymin);

freezenodes=abs(excrd(:,1)-xmin)<tolerance;
freezenodes=freezenodes | abs(excrd(:,1)-xmax)<tolerance;
freezenodes=freezenodes | abs(excrd(:,2)-ymin)<tolerance;
freezenodes=freezenodes | abs(excrd(:,2)-ymax)<tolerance;

% compute nodal connectivity
nodconnect=zeros(numnod,2);
for i=1:numseg
    nodeA = exedge(i, 1);
    nodeB = exedge(i, 2);
    if nodconnect(nodeA,1)==0; k=1; else; k=2; end
    nodconnect(nodeA,k) = i;
    if nodconnect(nodeB,1)==0; k=1; else; k=2; end
    nodconnect(nodeB,k) = i;
end

if pltflg>0
    figure(8)
    for i=1:numseg
        ni=exedge(i,1);
        nj=exedge(i,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);
        plot(xp,yp,'k.'); hold on;
    end
    for i=1:numnod
        if freezenodes(i)
            plot(excrd(i,1),excrd(i,2),'ro');hold on
        end
    end
    axis equal
end

inivec=ones(numnod,1);
inivec(freezenodes)=0;
dexcrd=diag(inivec);
newexcrd=excrd;
dnewexcrd=dexcrd;
for is=1:smonum
    for ip=1:2
        for in=1:numnod
            if freezenodes(in); continue; end

            n1=nodconnect(in,1);
            n2=nodconnect(in,2);

            newexcrd(in,:)=(1-smostep(ip))*excrd(in,:)+smostep(ip)/2*(excrd(n1,:)+excrd(n2,:));
            dnewexcrd(in,:)=(1-smostep(ip))*dexcrd(in,:)+smostep(ip)/2*(dexcrd(n1,:)+dexcrd(n2,:));
        end
        excrd=newexcrd;
        dexcrd=dnewexcrd;
    end
end

if pltflg>0
    figure(8)
    for i=1:numseg
        ni=exedge(i,1);
        nj=exedge(i,2);
        xp=excrd([ni nj],1);
        yp=excrd([ni nj],2);
        plot(xp,yp,'k-'); hold on;
    end
    axis equal
end
    
end