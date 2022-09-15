function [xyz,topo,ipa,fvf]=convertexo(exofname)
%
% extract data from exodus file
%
% input: exofname - filename of exodus file (string)
%
% output: xzy  - coordinates of nodes (numnodes x 3)
%         topo - connectivity of elements (numele x 8)
%         ipa  - inplane fiber angle (numele x 1)
%         fvf  - fiber volume fraction (numele x 1)
%
% example: [xyz,topo,ipa,fvf]=convertexo('cantbeam_matchd.exo');
%

% read exo file
data=ExtractExoData(0,exofname);

% determine number of nodes
numnodes=length(data.Coord{1}.Coords);

% determine number of element blocks
numblks=length(data.Etopo);

% determine number of elements
numele=0;
for ib=1:numblks
    numele=numele + size(data.Etopo{ib}.Nodes,2);
end

% initialize arrays for nodel coordiantes and element connectivity
xyz=zeros(numnodes,3);
topo=zeros(numele,8);

% extract nodal coordinates
for in=1:numnodes
    xyz(in,1)=data.Coord{1}.Coords(in);
    xyz(in,2)=data.Coord{2}.Coords(in);
    xyz(in,3)=data.Coord{3}.Coords(in);
end

% extract element connectivity
ice=0;
for ib=1:numblks
    for ie=1:size(data.Etopo{ib}.Nodes,2)
      ice=ice+1;
      topo(ice,:)=data.Etopo{ib}.Nodes(:,ie);
    end
end

% extract inplane angle: data.EVar{5,:}
ipa=zeros(numele,1);
ice=0;
for ib=1:numblks
    for ie=1:size(data.EVar{5,ib}.Val)
      ice=ice+1;
      ipa(ice)=data.EVar{5,ib}.Val(ie);
    end
end

% extract fiber volume fraction: data.EVar{6,:}
fvf=zeros(numele,1);
ice=0;
for ib=1:numblks
    for ie=1:size(data.EVar{6,ib}.Val)
      ice=ice+1;
      fvf(ice)=data.EVar{6,ib}.Val(ie);
    end
end
