function IniExoData(exofile,fecoord,thick)
% function IniExoData(exofile,fecoord)
% initializes exodus file 
%
% INPUT exofile : name of exodus file including extension in string format
%
%       fecoord : nodal coordinates 
%
% Sebastian Kreissl, University of Colorado, Boulder, Colorado
% 2010/Apr/02
%
% last modification 2010/Apr/02

% open exodus file
ncid = netcdf.open(exofile,'NC_WRITE');

% create pseudo 3D mesh coordinates
dimid=netcdf.inqDimID(ncid,'num_nodes');
[~, num_nodes] = netcdf.inqDim(ncid,dimid);

numfenod = size(fecoord,1);

% check if pseudo 3D mesh

if numfenod ~= num_nodes
    if num_nodes == 2*numfenod
        
        if nargin>2
            thk=thick;
        else
            thk=1; 
        end
        
        fectmp=zeros(2*numfenod,3);
        fectmp(1:numfenod,1:2)=fecoord;
        fectmp(numfenod+1:end,1:2)=fecoord;
        fectmp(numfenod+1:end,3)=thk;
        fecoord=fectmp;
    else
        error('inconsistent node numbers');
    end
end

% set coordinate names
switch size(fecoord,2)
    case 1
        CoordName{1} = 'coordx';
    case 2
        CoordName{1} = 'coordx';
        CoordName{2} = 'coordy';
    case 3
        CoordName{1} = 'coordx';
        CoordName{2} = 'coordy';
        CoordName{3} = 'coordz';
end


% switch to 'Def'-mode 
netcdf.reDef(ncid);

% change global attribute 'file_size' to 1
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'file_size',1);

%%%%%%%%%%%%%%%%%%%%
% define coordinates
%%%%%%%%%%%%%%%%%%%%
% exodus(dimension 4) is num_dim
[~, num_spacial_dims] = netcdf.inqDim(ncid,4);
spacial_dims_id = zeros(num_spacial_dims,1);
for isd = 1:1:num_spacial_dims
    % exodus(dimension 5) is num_nodes
    spacial_dims_id(isd) = netcdf.defVar(ncid,CoordName{isd},'double',5);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % define dimension for nodal variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:length(VarName)
%     str = ['vals_nod_var',num2str(i)];
%     % exodus(dimension 3) is time_step
%     netcdf.defVar(ncid,str,'double',[5,3]);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % define variable 'name_nod_var'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lvn = netcdf.defDim(ncid,'num_nod_var',length(VarName));
% name_nod_var_id = netcdf.defVar(ncid,'name_nod_var','char',[0,lvn]);


% leave 'Def'-mode
netcdf.endDef(ncid);

%%%%%%%%%%%%%%%%%%
% fill coordinates
%%%%%%%%%%%%%%%%%%
for isd = 1:1:num_spacial_dims
    netcdf.putVar(ncid,spacial_dims_id(isd),fecoord(:,isd));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write variable names to 'name_nod_var'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ci = 1:1:length(VarName)
%     for ri = 1:1:length(VarName{ci})
%         netcdf.putVar(ncid,name_nod_var_id,[ri-1,ci-1],VarName{ci}(ri));
%     end
% end

% close exodus file
netcdf.close(ncid)
