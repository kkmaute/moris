%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% REMEMBER TO DELETE LAST EMPTY LINE FROM DISTTABLE FILE!!!
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

% input

% Define variables, folderName without / at the end
folderName='/home/villanueva/work/optimization/structural/beam3D_structuralCompliance_hexa8_building_torsion_SIMP_massRatio_0.10_strainEnergy_126x030x030_nodesets_radius/';
meshName='beam3D_126x030x030.g';
rmax=0.16;  % Smoothing radius
numlev=3;   % Levels of adjacent elements to search for nodes in smoothing radius 
nproc=6;    % Number of processors to use in the search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add backslash to folder name
folderName=[folderName,'/'];

% Append folder name to mesh name
meshName=[folderName,meshName];

% Delete previous tables
cmd=['rm -f ',folderName,'table_*'];
system(cmd);

% Read mesh
[fecoord,fetopo,~,~,numfenod,numfele] = LoadExoMesh(meshName);

% Number of elems by three dimensions
centroid = zeros(numfele,3);

% Loop over all elements
for ielem = 1:numfele
    % Get number of nodes per element
    numNodesinElem = size(fetopo,2)-1;

    % Number of nodes by three dimensions
    nodeCoords = zeros(numNodesinElem,3);

    % Loop over the nodes in this element
    for inode = 1:numNodesinElem
        % Get node ID (first column is element type
        nodeID = fetopo(ielem,inode+1);
        
        % Get x,y,z coords
        nodeCoords(inode,:) = fecoord(nodeID,:);
    end
    
    % Find maximum and minimum to find centroid
    minX = min(nodeCoords(:,1));
    maxX = max(nodeCoords(:,1));
    minY = min(nodeCoords(:,2));
    maxY = max(nodeCoords(:,2));
    minZ = min(nodeCoords(:,3));
    maxZ = max(nodeCoords(:,3));
    
    % Centroid coordinates
    centroid(ielem,:) = [minX + (maxX-minX)/2 minY + (maxY-minY)/2 minZ + (maxZ-minZ)/2]; 
end

fname=sprintf([folderName,'disttable']);

fid=fopen(fname,'w');

% Loop over all elements in mesh
for ielem = 1:numfele
    
    % Checking progress
    fprintf('%f\n',ielem/numfele);
    if (ielem/numfele == 0.000250)
        fprintf('stop\n');
    end
    counter = 1;

    clear connectedNodeIDs;
    clear connectedNodeDstncs;

    % Loop over all nodes in mesh
    for inode = 1:numfenod
        nodeCoords = fecoord(inode,:);
        
        % Compute distance
        distance = sqrt((centroid(ielem,1) - nodeCoords(1,1))^2 + ...
                        (centroid(ielem,2) - nodeCoords(1,2))^2 + ...
                        (centroid(ielem,3) - nodeCoords(1,3))^2 );

        % See if inside smoothing radius
        if ( distance <= rmax )
            connectedNodeIDs(counter) = inode;
            connectedNodeDstncs(counter) = rmax - distance;
            
            counter = counter + 1;
        end
    end
    
    numNodesInRadius = counter - 1;
    assert(numNodesInRadius == size(connectedNodeIDs,2));
    
    normalization = sum(connectedNodeDstncs);
    
    % Print results to file
    fprintf(fid,'%1.7e %1.7e ',ielem,numNodesInRadius);
    for ii=1:numNodesInRadius
        fprintf(fid,'%1.7e ',connectedNodeIDs(ii));
    end
    for ii=1:numNodesInRadius
        fprintf(fid,'%1.7e ',connectedNodeDstncs(ii)/normalization);
    end
    
    if (ielem ~= numfele)
        fprintf(fid,'\n');
    end
end

