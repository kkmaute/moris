%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hardcoded for 2D flow!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
% clear all
% close all

format long e

% input

% Define variables, folderName without / at the end
folderName = '/home/villanueva/work/optimization/fefluid/turbine_Brinkman_research_smallMesh_symmetric_2time';
meshName = 'turbine.g';
rmax = 0.2;  % turbine_Brinkman_research_smallMesh
% rmax = 5.740089151084676e-02;  % turbine_Brinkman_research_smallMesh_finer

numlev = 3;   % Levels of adjacent elements to search for nodes in smoothing radius 
nproc = 6;    % Number of processors to use in the search

symmetry = true;

% Center of rotation
centerRotX = 0.0;
centerRotY = 0.0;

% Domain of rotation
domainRotXMin = -3.0;
domainRotXMax =  3.0;
domainRotYMin = -3.0;
domainRotYMax =  3.0;

% Rotational speed in radian per second
omega = 2 * pi;

% We want a quarter of a rotation
total_time = 1;
time_step = 0.01;
num_time_its = total_time / time_step;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read mesh file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add backslash to folder name
folderName=[folderName,'/'];

% Append folder name to mesh name
meshName=[folderName,meshName];

% Create folder if nonexistent
cmd=['mkdir  ',folderName,'disttable'];
system(cmd);

% Delete previous tables
cmd=['rm -f ',folderName,'disttable/table_*'];
system(cmd);

% Read mesh
[fecoord,fetopo,~,~,numfenod,numfele] = LoadExoMesh(meshName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find centroid of all elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of elems by three dimensions
centroid = zeros(numfele,2);

% Loop over all elements
parfor ielem = 1:numfele
    % Get number of nodes per element
    numNodesinElem = size(fetopo,2)-1;

    % Number of nodes by two dimensions
    nodeCoords = zeros(numNodesinElem,2);

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
%	minZ = min(nodeCoords(:,3));
%	maxZ = max(nodeCoords(:,3));

    % Centroid coordinates
% 	centroid(ielem,:) = [minX + (maxX-minX)/2 minY + (maxY-minY)/2 minZ + (maxZ-minZ)/2];
	centroid(ielem,:) = [minX + (maxX-minX)/2 minY + (maxY-minY)/2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it = 0:num_time_its

	tic
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Where was I before this?

    % Rotation
    time = time_step * it;
    theta = - omega * time;

    % Where was I before
    coordsPrevTime = zeros(numfele,2);
    newCentroid = zeros(numfele,2);

	% Compute where the elements were at this time step using vectors
	newCentroid(:,1) = centroid(:,1) - centerRotX;
	newCentroid(:,2) = centroid(:,2) - centerRotY;

	coordsPrevTime(:,1) =  cos(theta)*newCentroid(:,1) + sin(theta)*newCentroid(:,2);
	coordsPrevTime(:,2) = -sin(theta)*newCentroid(:,1) + cos(theta)*newCentroid(:,2);

	coordsPrevTime(:,1) = coordsPrevTime(:,1) + centerRotX;
	coordsPrevTime(:,2) = coordsPrevTime(:,2) + centerRotY;
	
	% Model ADVs symmetrically
	if ( symmetry == true )
		for ielem = 1:numfele
			
			% Domain 1 do nothing

			% Domain 2: rotate CCW 90 degrees
			if     ( coordsPrevTime(ielem,1) >  0.0 && coordsPrevTime(ielem,2) >= 0.0 )
				newCentroid(ielem,1) = coordsPrevTime(ielem,1) - centerRotX;
				newCentroid(ielem,2) = coordsPrevTime(ielem,2) - centerRotY;
				
				coordsPrevTime(ielem,1) =  cosd(- 90)*newCentroid(ielem,1) + sind(- 90)*newCentroid(ielem,2);
				coordsPrevTime(ielem,2) = -sind(- 90)*newCentroid(ielem,1) + cosd(- 90)*newCentroid(ielem,2);
				
				coordsPrevTime(ielem,1) = coordsPrevTime(ielem,1) + centerRotX;
				coordsPrevTime(ielem,2) = coordsPrevTime(ielem,2) + centerRotY;
				
% 				coordsPrevTime(ielem,1) = -1.0 * coordsPrevTime(ielem,1);

			% Domain 3: rotate CCW 180 degrees
			elseif ( coordsPrevTime(ielem,1) >  0.0 && coordsPrevTime(ielem,2) <  0.0 )

				coordsPrevTime(ielem,1) = -1.0 * coordsPrevTime(ielem,1);
				coordsPrevTime(ielem,2) = -1.0 * coordsPrevTime(ielem,2);

			% Domain 4: rotate CCW 270 degrees
			elseif ( coordsPrevTime(ielem,1) <= 0.0 && coordsPrevTime(ielem,2) <  0.0 )
				newCentroid(ielem,1) = coordsPrevTime(ielem,1) - centerRotX;
				newCentroid(ielem,2) = coordsPrevTime(ielem,2) - centerRotY;
				
				coordsPrevTime(ielem,1) =  cosd(-270)*newCentroid(ielem,1) + sind(-270)*newCentroid(ielem,2);
				coordsPrevTime(ielem,2) = -sind(-270)*newCentroid(ielem,1) + cosd(-270)*newCentroid(ielem,2);
				
				coordsPrevTime(ielem,1) = coordsPrevTime(ielem,1) + centerRotX;
				coordsPrevTime(ielem,2) = coordsPrevTime(ielem,2) + centerRotY;

% 				coordsPrevTime(ielem,2) = -1.0 * coordsPrevTime(ielem,2);
			end
		end
	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute weights

    fname=sprintf([folderName,'disttable/disttable_','%1.3f'],time);

    fid=fopen(fname,'w');

    % Loop over all elements in mesh
    for ielem = 1:numfele

% 		if ielem == 685
% 			fprintf('Stop\n');
% 		end

%		counter = 1;

        clear connectedElemIDs;
        clear connectedElemIDs2;
        clear connectedElemDstncs;

        connectedElemIDs = [];
        connectedElemIDs2 = [];
        connectedElemDstncs = [];
		
		% If inside element of rotation
% 		if ( centroid(ielem,1) >= domainRotXMin && centroid(ielem,1) <= domainRotXMax && centroid(ielem,2) >= domainRotYMin && centroid(ielem,2) <= domainRotYMax )
			
        % Loop over all nodes in mesh
%       for jelem = 1:numfele

            % Compute distance
    %       distance = sqrt((centroid(ielem,1) - coordsPrevTime(jelem,1))^2 + ...
    %                       (centroid(ielem,2) - coordsPrevTime(jelem,2))^2 + ...
    %                       (centroid(ielem,3) - coordsPrevTime(jelem,3))^2 );
            distance = sqrt((centroid(:,1) - coordsPrevTime(ielem,1)).^2 + ...
                            (centroid(:,2) - coordsPrevTime(ielem,2)).^2);

            connectedElemIDs = find(distance <= rmax);
			
			if ( symmetry )
				counter = 1;
				for m = 1:size(connectedElemIDs,1)
					if ( centroid(connectedElemIDs(m),1) <= 0.0 && centroid(connectedElemIDs(m),2) >= 0.0 && ...
							sqrt( ( centroid(connectedElemIDs(m),1) - centerRotX ) ^ 2 + ( centroid(connectedElemIDs(m),2) - centerRotY ) ^ 2 ) <= 3.0 )
						connectedElemIDs2(counter,1) = connectedElemIDs(m);
						counter = counter + 1;
					end
				end
			end
			connectedElemIDs = connectedElemIDs2;
            connectedElemDstncs = rmax - distance(connectedElemIDs,1);

            % See if inside smoothing radius
%         	if ( distance <= rmax )
%           	connectedElemIDs(counter) = jelem;
%               connectedElemDstncs(counter) = rmax - distance;
% 
%           	counter = counter + 1;
%       	end
%       end

% 		end

        numElemsInRadius = size(connectedElemIDs,1);

%   	numElemsInRadius = counter - 1;
%     	assert(numElemsInRadius == size(connectedElemIDs,2));

		normalization = sum(connectedElemDstncs);

        % Print results to file
        fprintf(fid,'%1.7e %1.7e ',ielem,numElemsInRadius);
        for ii=1:numElemsInRadius
            fprintf(fid,'%1.7e ',connectedElemIDs(ii));
        end
        for ii=1:numElemsInRadius
            fprintf(fid,'%1.7e ',connectedElemDstncs(ii)/normalization);
        end

        if (ielem ~= numfele)
            fprintf(fid,'\n');
        end
	end
	
	% Checking progress
	fprintf('%f\n',it/num_time_its);

	toc
end
