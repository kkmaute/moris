clear
close all
clc
format long e

%% Input.

scripts_dir = 'C:\Users\maute\Documents\MATLAB\TransFem Current\Matlab\misc';
problem_dir = 'C:\Users\maute\Desktop';
datafile  = 'fuzzy.e-s.0125';
meshfile  = 'fuzzy_refined.g'
rhobound  = 0.01;
refstrat  = 1;
iblk      = 2;

%% Read data.

owd = pwd;
addpath(scripts_dir);

%% Setup.

cd(problem_dir);

%% Extract information.

Data = ExtractExoData(0, datafile);

if size(Data.Coord,2)
    xdim=2;
    fecoord = [Data.Coord{1}.Coords Data.Coord{2}.Coords];
else
    xdim=3;
    fecoord = [Data.Coord{1}.Coords Data.Coord{2}.Coords Data.Coord{3}.Coords];
end

%% Process only element block defined by iblk

topo = Data.Etopo{iblk}.Nodes;

for i = 1 : size(Data.EVar, 1)
    if strfind(Data.EVar{i,iblk}.Name, 'MatPropStrucDensity') == 1
        density_evar = i;
        break;
    end
end

erho = Data.EVar{density_evar,iblk}.Val(:, 1);

numfele = size(topo, 2);

%% Compute elemental and nodal level set values for blockset 1.

switch refstrat
    
    case 1 % based on elemental density value
        
        tobeRemoved = find(erho<rhobound);
        
    otherwise
        error('incorrect refinement strategy');
end

%% Create refinement script.

fid = fopen('refine_simp.txt', 'w');
fprintf(fid, 'reset\n');
fprintf(fid, 'import mesh geometry "%s" block %d use nodeset sideset feature_angle 135.00 linear merge\n', [problem_dir, '/', datafile],iblk);
fprintf(fid, 'set dev on\n');

fprintf('number of elements to be deleted: %d\n\n', length(tobeRemoved));

%     if (xdim==2)
%         fprintf(fid, 'set duplicate block elements on\nblock 111 face ');
%     else
%         fprintf(fid, 'set duplicate block elements on\nblock 111 hex ');
%     end
%
%     for ir = 1 : length(tobeRefined)
%         fprintf(fid, '%d,', tobeRefined(ir));
%     end

if (xdim==2)
    fprintf(fid, 'delete face ');
else
    fprintf(fid, 'refine hex ');
end

for ir = 1 : length(tobeRemoved)
    fprintf(fid, '%d,', tobeRemoved(ir));
end

if (xdim==2)
    fprintf(fid, '\nrefine surface all ');
else
    fprintf(fid, '\nrefine volume all ');
end

%     for ir = 1 : length(tobeRefined)
%         fprintf(fid, '%d,', tobeRefined(ir));
%     end
%
fprintf(fid, 'numsplit 1 bias 1 depth 1 \n');
fprintf(fid, 'export mesh "%s" dimension %i overwrite\n', meshfile, xdim);

fclose(fid);

