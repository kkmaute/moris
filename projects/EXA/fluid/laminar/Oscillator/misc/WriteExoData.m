function WriteExoData(exofile,NodVarName,NodData,TimeStep,Time,EleVarName,EleData)
% function WriteExoData(exofile,NodVarName,NodData,TimeStep,Time)
% writes nodal (optional elemental) data and and time to exodus file
%
% INPUT exofile: name of exodus file including extension in string format
%
%       NodVarName  : names of variables in cell strings
%                  eg: NodVarName{1} = 'ux';
%                      NodVarName{2} = 'uy;
%                      NodVarName{3} = 'pres'
%
%       NodData: nodal data in colomn format
%                 eg: [ux_node1, uy_node1, pres_node1;
%                      ux_node2, uy_node2, pres_node2];
%
%       TimeStep: current time-step (INTEGER)
%
%       Time    : current time (REAL)
%
%       EleVarName  : names of variables in cell strings
%                  eg: EleVarName{1} = 'Porosity';
%                      EleVarName{2} = 'YoungsMod';
%
%       EleData: elemental data in colomn format
%                 eg: [Porosity_ele1, YoungsMod_ele1;
%                      Porosity_ele2, YoungsMod_ele2];
%
% Sebastian Kreissl, University of Colorado, Boulder, Colorado
% 2010/Apr/02
%
% last modification 2011/Feb/09

[LN,sND] = size(NodData);
lVN = length(NodVarName);

% open exodus file
ncid = netcdf.open(exofile,'NC_WRITE');

% get number of dimensions and number of variables from exodus file
[~,nvars] = netcdf.inq(ncid);

dimid = netcdf.inqDimID(ncid,'num_el_blk');
[~, num_blocks] = netcdf.inqDim(ncid,dimid);

time_step_id=netcdf.inqDimID(ncid,'time_step');

num_nodes_id=netcdf.inqDimID(ncid,'num_nodes');

[~, num_nodes] = netcdf.inqDim(ncid,num_nodes_id);

% check if pseudo 3D mesh
if LN ~= num_nodes
    if num_nodes == 2*LN
        NodData = [NodData;NodData];
        LN      = 2*LN;
    else
        error('inconsistent node numbers');
    end
end

% check if number of variable names matches with number of variable data
if sND~=lVN
    error('number of nodal variable names does not match with number of nodal variable data');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if variable names exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nod_var_exist_flag = zeros(lVN,1);
num_nod_var_exist_flag = 0;
for varid = 0:1:nvars-1
    varname = netcdf.inqVar(ncid,varid);
    for i = 1:1:lVN
        if strcmp(varname,'num_nod_var')
            num_nod_var_exist_flag = 0;
        end
        str = ['vals_nod_var',num2str(i)];
        if strcmp(varname,str)
            nod_var_exist_flag(i) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create nodal variables that do not exist yet
% => this should happen only the first time this function is called
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(nod_var_exist_flag) ~= sND
    
    % switch to 'Def'-mode
    netcdf.reDef(ncid);
    
    if num_nod_var_exist_flag == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define variable 'name_nod_var'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lvn = netcdf.defDim(ncid,'num_nod_var',lVN);
        name_nod_var_id = netcdf.defVar(ncid,'name_nod_var','char',[0,lvn]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define dimension for nodal variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:1:lVN
        if nod_var_exist_flag(i) == 0
            
            str = ['vals_nod_var',num2str(i)];
            netcdf.defVar(ncid,str,'double',[num_nodes_id,time_step_id]);
            
        end
    end
    
    % leave 'Def'-mode
    netcdf.endDef(ncid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write variable names to 'name_nod_var'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:1:lVN
        if nod_var_exist_flag(i) == 0
            
            for ci = 1:1:length(NodVarName)
                for ri = 1:1:length(NodVarName{ci})
                    netcdf.putVar(ncid,name_nod_var_id,[ri-1,ci-1],NodVarName{ci}(ri));
                end
            end
            
        end
        
    end
    
end

%%%%%%%%%%%
% fill Time
%%%%%%%%%%%
varid = netcdf.inqVarID(ncid,'time_whole');
netcdf.putVar(ncid,varid,TimeStep-1,Time);

%%%%%%%%%%%%%%%%%%%%%%
% fill nodal variables
%%%%%%%%%%%%%%%%%%%%%%
for ci = 1:1:sND
    % get ID of variable by name
    str = ['vals_nod_var',num2str(ci)];
    varid = netcdf.inqVarID(ncid,str);
    netcdf.putVar(ncid,varid,[0,TimeStep-1],[LN,1],NodData(:,ci));
end

%%%%%%%%%%%%%%%%
% elemental data
%%%%%%%%%%%%%%%%
if nargin > 5
    
    [LE,sED] = size(EleData);
    lVE = length(EleVarName);
    
    ele_per_block = zeros(num_blocks,1);
    ele_per_block_dimid = zeros(num_blocks,1);
    for i = 1:1:num_blocks
        str = ['num_el_in_blk',num2str(i)];
        ele_per_block_dimid(i) = netcdf.inqDimID(ncid,str);
        [~, ele_per_block(i)] = netcdf.inqDim(ncid,ele_per_block_dimid(i));
        if i == 1
            list = 1:ele_per_block(i);
        else
            list =  ele_per_block(i-1)+1:ele_per_block(i-1)+ele_per_block(i);
        end
        EleData_vec{i} =  EleData(list,:);
    end
    
    % check if number of variable names matches with number of variable data
    if sED~=lVE
        error('number of elemental variable names does not match with number of elemental variable data');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check if variable names exist
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elem_var_exist_flag = zeros(lVE,1);
    num_elem_var_exist_flag = 0;
    for varid = 0:1:nvars-1
        varname = netcdf.inqVar(ncid,varid);
        for i = 1:1:lVE
            if strcmp(varname,'num_elem_var')
                num_elem_var_exist_flag = 0;
            end
            str = ['vals_elem_var',num2str(i),'eb1'];   % check only for block 1
            if strcmp(varname,str)
                elem_var_exist_flag(i) = 1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create elemental variables that do not exist yet
    % => this should happen only the first time this function is called
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(elem_var_exist_flag) ~= sED
        
        % switch to 'Def'-mode
        netcdf.reDef(ncid);
        
        if num_elem_var_exist_flag == 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % define variable 'name_eleme_var'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lve = netcdf.defDim(ncid,'num_elem_var',lVE);
            name_elem_var_id = netcdf.defVar(ncid,'name_elem_var','char',[0,lve]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define dimension for elemental variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:1:num_blocks
            for i = 1:1:lVE
                if elem_var_exist_flag(i) == 0
                    
                    str = ['vals_elem_var',num2str(i),'eb',num2str(j)];
                    % exodus(dimension 3) is time_step
                    netcdf.defVar(ncid,str,'double',[ele_per_block_dimid(j),3]);
                    
                end
            end
        end
        
        % leave 'Def'-mode
        netcdf.endDef(ncid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write variable names to 'name_elem_var'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:1:lVE
            if elem_var_exist_flag(i) == 0
                
                for ci = 1:1:length(EleVarName)
                    for ri = 1:1:length(EleVarName{ci})
                        netcdf.putVar(ncid,name_elem_var_id,[ri-1,ci-1],EleVarName{ci}(ri));
                    end
                end
                
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill elemental variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for bi = 1:1:num_blocks
        for ci = 1:1:sED
            % get ID of variable by name
            str = ['vals_elem_var',num2str(ci),'eb',num2str(bi)];
            varid = netcdf.inqVarID(ncid,str);
            netcdf.putVar(ncid,varid,[0,TimeStep-1],[ele_per_block(bi),1],EleData_vec{bi}(:,ci));
        end
    end
end

% close exodus file
netcdf.close(ncid)

