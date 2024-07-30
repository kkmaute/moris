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

% flag for ploting structure
pltflg=1;

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
ic=0;
for i=1:numquads
    nids=[quads(i,:) quads(i,1)];
    for j=1:4
        ni=quads(i,j);
        nj=nids(j+1);
        ic=ic+1;
        edgcon(ic,:)=[sort([ni nj],'ascend') i];
    end
end
for i=1:numtris
    nids=[tris(i,:) tris(i,1)];
    for j=1:3
        ni=tris(i,j);
        nj=nids(j+1);
        ic=ic+1;
        edgcon(ic,:)=[sort([ni nj],'ascend') i+numquads];
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