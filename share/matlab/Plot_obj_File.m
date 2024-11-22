clear
clc
close all

%% User Inputs

% Surface mesh inputs
fileName = "";   % string for file name including full path but NOT extension
iterations = 0;                                                 % Optimization iteration you wish to plot

% FEM Mesh inputs
numElems = [2;2;2];                 % Vector of length dims containing the number of elements for each spatial dimension
domainDims = [-1 1; -1 1; -1 1;];   % dims x 2 matrix where the first column is lower bound and the second column is upper bound

%% Obj File Reader
figure();
hold on;
grid on;
for i = iterations
    % Get the file name and open it for reading
    file = strcat(fileName, '_', num2str(i), '.txt');
    FID = fopen(file, 'r');

    % initialize counters for number of vertices and facets
    numVertices = 0;
    numFacets = 0;
    
    % Get the first line of the file and split it into its components
    line = fgetl(FID);
    firstLine = split(line, ' ');
    
    % Ensure there are no empty elements at the end of the string
    for j = length(firstLine):-1:0
        if ~isempty(cell2mat(firstLine(j)))
            dims = j-1;
            break;
        end
    end

    % Loop through rest of the file and count vertices and facets
    while ischar(line)
        if( ~isempty(line) )
            if( strcmp(line(1:2), 'v ') )
                numVertices = numVertices + 1;
    
            elseif( strcmp(line(1:2), 'f ') )
                numFacets = numFacets + 1;
            end
        end

        line = fgetl(FID);
    end

    % Initialize vertex and facet arrays
    vertices = zeros(numVertices, dims);
    facets = zeros(numFacets, dims);

    % Reset counters
    numVertices = 1;
    numFacets = 1;

    % Loop through the file and read the vertices and facets
    frewind(FID);
    line = fgetl(FID);
    while ischar(line)
        % Get the parts of the string separated by spaces
        strElems = split(line, ' ');
        if ( strElems{1} == 'v' )
            vertices(numVertices, :) = cellfun(@str2num, strElems(2:dims+1));
            numVertices = numVertices + 1;

        elseif ( strElems{1} == 'f' )
            facets(numFacets, :) = cellfun(@str2num, cellfun(@parseFacetLine, strElems(2:dims+1), 'UniformOutput', false ));
            numFacets = numFacets + 1;
        end

        line = fgetl(FID);
    end

    % Plot this surface mesh
    if dims == 3
        trisurf(facets, vertices(:,1), vertices(:,2), vertices(:,3), 'FaceAlpha', .9 );
    elseif dims == 2
        for j=1:size(facets, 1)
            x1 = vertices(facets(j, 1), 1);
            x2 = vertices(facets(j, 2), 1);
            y1 = vertices(facets(j, 1), 2);
            y2 = vertices(facets(j, 2), 2);
            plot([x1 x2], [y1 y2])
            % fprintf('x: %f\ty:%f\n', x(i), y(i));
        end
    end

    fclose(FID);
end

% Set the FEM Mesh up as the grid
xlim(domainDims(1,:));
ylim(domainDims(2,:));

xticks(linspace(domainDims(1,1), domainDims(1,2), numElems(1) + 1 ));
yticks(linspace(domainDims(2,1), domainDims(2,2), numElems(2) + 1 ));

xlabel('x');
ylabel('y');

if dims == 3
    zlim(domainDims(3,:));
    zticks(linspace(domainDims(3,1), domainDims(3,2), numElems(3) + 1 ));
    zlabel('z');
end

%% Function definitions
function parsed = parseFacetLine(s)
    if ~contains(s, '/')
        parsed = s;
        return;

    else 
        indicesOfDelimiter = strfind(s, '/');
        parsed = s(1:indicesOfDelimiter(1)-1);
    end
end