function [fecoord,fetopo,feprop,febc,numfenod,numfele]= LoadExoMesh(meshname,numblkused)
% function [fecoord,fetopo,feprop,febc,numfenod,numfele]= LoadExoMesh(meshname)
% reads a cubit meshfile and creates /initializes and exodus file for
% paraview
% 
% INPUT: meshname: name of mesh in string format with file extension
%                  e.g. LoadExoMesh('ExampleMesh.g')
%
% OUTPUT:
%       fecoord   - coordinates of fe nodes (numfenod x 2)
%       fetopo    - element topology (4-node quads) (numfele x 5)
%                   fetopo(1,:) : 1 = volume element
%                                 2 = boundary element (belongs to side set)
%       feprop    - element properties (numfele x 2)
%                   1. column: 0 = volume element
%                              0 < side set index
%                   2. column  block set id
%       febc      - node set or side set id (numfenod x 1)
%                   watch: nodes with conflicting node/side set ids get
%                          information from last node/side set 
%                          (not necesseraly id)
%                   watch: side set id overwrites node set id
%       numfenod    number of fe nodes
%       numfele     number of elements
%

% flags for element types in fetopo
vol_ele_flag = 1;
bndry_ele_flag = 2;

% print information from ExtractExoData: 0: no | 1: yes
printflag = 0;

% extract mesh information from exodus file
MeshData = ExtractExoData(printflag,meshname);

% assign mesh information

% nodal coordiantes 
for i=1:1:length(MeshData.Coord)
    fecoord(:,i) = MeshData.Coord{1,i}.Coords;
end

% elements stored in blocks
numblk=size(MeshData.Etopo,2);
if numblkused < numblk
    numblk=numblkused;
end

% initialize elemental topology and properties
fetopo=[];
feprop=[];
feid=[];

numfele=0;

for ib=1:numblk

    fetopo_tmp  = MeshData.Etopo{1,ib}.Nodes';
    numfele_tmp = size(fetopo_tmp,1);

    fetopo = [fetopo;vol_ele_flag*ones(numfele_tmp,1),fetopo_tmp];
    feprop = [feprop;zeros(numfele_tmp,1),ib*ones(numfele_tmp,1)];
    
    if isempty(MeshData.EleMap)
        elemap=numfele+1:numfele+numfele_tmp;
    else
        elemap=MeshData.EleMap(numfele+1:numfele+numfele_tmp);
    end
    
    feid   = [feid,elemap'];
    numfele = numfele+numfele_tmp;
end

% resort elements
fetopo=fetopo(feid,:);
feprop=feprop(feid,:);

% total number of nodes
numfenod = size(fecoord,1);

% set boundary conditions
febc=zeros(numfenod,1);

% Side to Node map for Quad 4 element;
switch size(fetopo,2)
    case 4
        SideNodeMap = [1 2; 2 3; 3 1];
    case 5
        SideNodeMap = [1 2; 2 3; 3 4; 4 1];
    case 9
        SideNodeMap = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    otherwise
        error('element type not implemented');
end
        
% for given node sets
if isempty(MeshData.NodeSet) == 0

    febc_NS = zeros(numfenod,1);
    num_BC_NS = length(MeshData.NodeSet);
    
    for bi = 1:1:num_BC_NS
        febc_NS(MeshData.NodeSet{1,bi}.Nodes) = MeshData.NodeSetID{1}(bi);
    end
    
    % write the BC entries from febc_SS into febc
    
    ind_NS = find(febc_NS~=0);
    febc(ind_NS) = febc_NS(ind_NS);
end

% for given side sets
if isempty(MeshData.SideSet) == 0

    febc_SS = zeros(numfenod,1);
    num_BC_SS = length(MeshData.SideSet);
    
    for bi = 1:1:num_BC_SS

        ele_per_SS = length(MeshData.SideSet{1,bi}.Elements);

        for ei = 1:1:ele_per_SS

            % flag boundary elements with "2" in fetopo
            ele_ind = MeshData.SideSet{1,bi}.Elements(ei);
            fetopo(ele_ind,1) = bndry_ele_flag;
            
            % assign the side index of the boundary element to feprop
            side_ind = MeshData.SideSet{1,bi}.Sides(ei);
            feprop(ele_ind,1) = side_ind;
            febc_SS(fetopo(ele_ind,SideNodeMap(side_ind,:))) =  MeshData.SideSetID{1}(bi);
        end
    end
    
    % write the BC entries from febc_SS into febc
    ind_SS = find(febc_SS~=0);
    febc(ind_SS) = febc_SS(ind_SS);
end