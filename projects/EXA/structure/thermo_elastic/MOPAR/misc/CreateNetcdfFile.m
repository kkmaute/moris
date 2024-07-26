function CreateNetcdfFile(filename,eletype,coordinate_names,fecoord,fetopo,numfele,febc,num_el_blk,num_el_in_blk,num_side_sets,sideset_ind,num_side_ss,num_df_ss,dist_fact_ss,ele_in_sideset,side_in_sideset)
% function % IniExoFile(filename,eletype,coordinate_names,fecoord,fetopo,numfele,febc,num_el_blk,num_el_in_blk,num_side_sets,sideset_ind,num_side_ss,num_df_ss,dist_fact_ss,ele_in_sideset,side_in_sideset)
% Creates an Netcdf file(*.g) from MATLAB data
%
% % Sebastian Kreissl, University of Colorado, Boulder, Colorado
% 2010/Apr/27
%
% last modification 2010/Apr/29
%
% INPUT filename :        name of netcdf file including extension
%       eletype:          a cell of strings one for every blockset
%                         e.g. eletype = {'quad4','quad8'}
%       coordinate_names: cell of coordinates eg. {'x','z'}
%       fecoord:          nodal coordinates, rows correspond to node,
%                         columns to spatial dimension
%       fetopo:           connectivity matrix nodes are
%                         dimension for quad4: (# of elements) x 5
%       numfele:          number of elements
%       febc:             column vector of boundary flags
%                         dimension (# of nodes) x 1
%       num_el_blk:       number of blocks
%                         dimension: (# of blocks) x 1
%       num_el_in_blk:    number of elements in blocks
%                         dimension: (# of blocks) x 1
%       num_side_sets:    number of sidesets
%       sideset_ind:      (optional) indices of sidesets
%                         dimension: (# of sidesets) x 1
%       num_side_ss:      (optional) number of sides per sideset
%                         dimension: (# of sideset,1)
%                         e.g.: num_side_ss = [2,2] => 2 sidessets with two
%                         sides each
%       num_df_ss:        (optional) number of distribution factors (each
%                         side has two (?) => 2*num_side_ss;
%       dist_fact_ss:     (optional): right now all are assumed to be one
%                         => dist_fact_ss = ones(# number of sidesets,2)
%       ele_in_sideset:   (optional) elements in side sets
%                         e.g.: ele_in_sideset = [1 2 0; 4 5 6]
%                         => sideset #1 contains elements 1, 2
%                            sideset #2 contains elements 4, 5, 6
%       side_in_sideset:  (optional) elements in side sets
%                         analogous to ele_in_sideset
%
%

% create quasi 3D mesh
if strfind(filename,'_3D')
    if strcmp(eletype,'quad4')
        numfenod = size(fecoord,1);
        numfele=size(fetopo,1);
        fectmp=zeros(2*numfenod,3);
        ftptmp=zeros(numfele,9);
        
        fectmp(1:numfenod,1:2)=fecoord;
        fectmp(numfenod+1:end,1:2)=fecoord;
        fectmp(numfenod+1:end,3)=1;            % thickness
        ftptmp(:,2:5)=fetopo(:,2:5);
        ftptmp(:,6:9)=fetopo(:,2:5)+numfenod*ones(numfele,4);
        
        eletype={'hex8'};
        coordinate_names={'x','y','z'};
        fecoord=fectmp;
        fetopo=ftptmp;
    end
end
        
% number of nodes
numfenod = size(fecoord,1);
% number of spatial dimensions
num_dim = size(fecoord,2);

% determine node set flags
node_set_ind = unique(febc);
node_set_ind(node_set_ind == 0) = [];
ns_prop1 = node_set_ind;
num_node_sets = length(ns_prop1);

% num_nod_ns: how many nodes are in each node set
num_nod_ns = zeros(1,num_node_sets);
% num_df_ns: distribution factors
num_df_ns = zeros(1,num_node_sets);
for i = 1:1:num_node_sets
    num_nod_ns(i) = length(find(febc==node_set_ind(i)));
    num_df_ns(i) = num_nod_ns(i);
end

% distribution factors for nodes in nodesets (is not used now)
dist_fact_ns = ones(num_node_sets,max(num_nod_ns));
% no idea what that means => maybe to indicate non-empty nodesets
ns_status = ones(num_node_sets,1);

% nodes_in_node_set: row corresponds to node set, columns contain nodes
nodes_in_nodeset = zeros(num_node_sets,max(num_nod_ns));
for i = 1:1:num_node_sets
    nodes_in_nodeset(i,1:num_nod_ns(i)) = find(febc==node_set_ind(i));
end

% ATTENTION: no idea what these variables mean
eb_status_val = ones(1,num_el_blk); % no idea
eb_prop1_val = 1:num_el_blk; % no idea

if nargin <= 9
    % set number of sidesets to zero
    num_side_sets = 0;
else
    % no idea what that means => maybe to indicate non-empty sidesets
    ss_status = ones(num_side_sets,1);
    % side set flags
    ss_prop1 = sideset_ind;
end

title = 'none';

num_nod_per_el=zeros(num_el_blk,1);

for i = 1:1:num_el_blk
    switch eletype{i}
        case 'bar2'
            num_nod_per_el(i) = 2;
        case 'tri3'
            num_nod_per_el(i) = 3;
        case 'tri6'
            num_nod_per_el(i) = 6;
        case 'quad4'
            num_nod_per_el(i) = 4;
        case 'tet4'
            num_nod_per_el(i) = 4;            
        case 'tet10'
            num_nod_per_el(i) = 10;            
        case 'quad8'
            num_nod_per_el(i) = 8;
       case 'quad9'
            num_nod_per_el(i) = 9;
        case 'hex8'
            num_nod_per_el(i) = 8;
        case 'hex20'
            num_nod_per_el(i) = 20;            
        otherwise
            error('incorrect element type');
    end
end

int_str = 'int';
double_str = 'double';
char_str = 'char';

% create file
ncid = netcdf.create(filename, 'NC_SHARE');
netcdf.close(ncid);

% open file
ncid = netcdf.open(filename,'NC_WRITE');

% switch to 'Def'-mode
netcdf.reDef(ncid);

% create dimensions
len_string_dimID = netcdf.defDim(ncid,'len_string',33);
len_line_dimID = netcdf.defDim(ncid,'len_line',81);
four_dimID = netcdf.defDim(ncid,'four',4);
time_step_dimID = netcdf.defDim(ncid,'time_step',netcdf.getConstant('NC_UNLIMITED'));
num_dim_dimID = netcdf.defDim(ncid,'num_dim',num_dim);
num_nodes_dimID = netcdf.defDim(ncid,'num_nodes',numfenod);
num_elem_dimID = netcdf.defDim(ncid,'num_elem',numfele);
num_el_blk_dimID = netcdf.defDim(ncid,'num_el_blk',num_el_blk);
num_qa_rec_dimID = netcdf.defDim(ncid,'num_qa_rec',1);

% if node sets exist
if num_node_sets > 0
    num_node_sets_dimID = netcdf.defDim(ncid,'num_node_sets',num_node_sets);
    
    for i=1:1:num_node_sets
        % create num_nod_ns1, num_nod_ns2, etc.
        num_nod_ns_str = ['num_nod_ns',num2str(i)];
        name_num_nod_ns_dimID = genvarname(['num_nod_ns',num2str(i),'_dimID']);
        eval([name_num_nod_ns_dimID '= netcdf.defDim(ncid, num_nod_ns_str,num_nod_ns(i));']);
    end
    
end

% if side sets exist
if num_side_sets > 0
    num_side_sets_dimID = netcdf.defDim(ncid,'num_side_sets',num_side_sets);
    
    for i=1:1:num_side_sets
        % create num_side_ss1, num_side_ss2, etc.
        num_side_ss_str = ['num_side_ss',num2str(i)];
        name_num_side_ss_dimID = genvarname(['num_side_ss',num2str(i),'_dimID']);
        eval([name_num_side_ss_dimID '= netcdf.defDim(ncid, num_side_ss_str,num_side_ss(i));']);
        
        % create num_df_ss1, num_df_ss2, etc.
        num_df_ss_str = ['num_df_ss',num2str(i)];
        name_num_df_ss_str = genvarname(num_df_ss_str);
        eval([name_num_df_ss_str '=netcdf.defDim(ncid, num_df_ss_str,num_df_ss(i));']);
    end
end

for i =1:1:num_el_blk
    % create num_el_in_blk1, num_el_in_blk2, etc.
    num_el_in_blk_str = ['num_el_in_blk',num2str(i)];
    name_num_el_in_blk_dimID = genvarname(['num_el_in_blk',num2str(i),'_dimID']);
    eval([name_num_el_in_blk_dimID '= netcdf.defDim(ncid, num_el_in_blk_str,num_el_in_blk(i));']);
    
    % create num_node_per_el1, num_node_per_el2, etc
    num_nod_per_el_str = ['num_nod_per_el',num2str(i)];
    name_num_nod_per_el_dimID = genvarname(['num_nod_per_el',num2str(i),'_dimID']);
    eval([name_num_nod_per_el_dimID '= netcdf.defDim(ncid, num_nod_per_el_str,num_nod_per_el(i));']);
end

% create variables
time_whole_varID = netcdf.defVar(ncid,'time_whole','double',time_step_dimID);
qa_records_varID = netcdf.defVar(ncid,'qa_records','char',[len_string_dimID,four_dimID,num_qa_rec_dimID]);
coor_names_varID = netcdf.defVar(ncid,'coor_names','char',[len_string_dimID,num_dim_dimID]);
eb_names_varID = netcdf.defVar(ncid,'eb_names','char',[len_string_dimID,num_el_blk_dimID]);

if num_node_sets > 0
    ns_status_varID = netcdf.defVar(ncid,'ns_status','int',num_node_sets_dimID);
    ns_prop1_varID = netcdf.defVar(ncid,'ns_prop1','int',num_node_sets_dimID);
    netcdf.putAtt(ncid,ns_prop1_varID,'name','ID');
    ns_names_varID = netcdf.defVar(ncid,'ns_names','char',[len_string_dimID,num_node_sets_dimID]);
    
    for i=1:1:num_node_sets
        % set node_ns1, node_ns2, etc
        node_ns_str = ['node_ns',num2str(i)];
        name_node_ns_varID = genvarname(['node_ns',num2str(i),'_dimID']);
        num_nod_ns_str = ['num_nod_ns',num2str(i),'_dimID'];
        eval([name_node_ns_varID '= netcdf.defVar(ncid,node_ns_str,int_str,eval(num_nod_ns_str));']);
        
        % set dist_fact_ns1, dist_fact_ns2, etc
        dist_fact_ns_str = ['dist_fact_ns',num2str(i)];
        num_nod_ns_str = ['num_nod_ns',num2str(i),'_dimID'];
        name_dist_fact_ns_varID = genvarname(['dist_fact_ns',num2str(i),'_dimID']);
        eval([name_dist_fact_ns_varID '= netcdf.defVar(ncid,dist_fact_ns_str,double_str,eval(num_nod_ns_str));']);
    end
end

if num_side_sets > 0
    ss_status_varID = netcdf.defVar(ncid,'ss_status','int',num_side_sets_dimID);
    ss_prop1_varID = netcdf.defVar(ncid,'ss_prop1','int',num_side_sets_dimID);
    netcdf.putAtt(ncid,ss_prop1_varID,'name','ID');
    ss_names_varID = netcdf.defVar(ncid,'ss_names','char',[len_string_dimID,num_side_sets_dimID]);
    
    for i=1:1:num_side_sets
        % set elem_ss1, elem_ss2, etc
        elem_ss_str = ['elem_ss',num2str(i)];
        name_elem_ss_varID = genvarname(['elem_ss',num2str(i),'_dimID']);
        num_side_ss_str = ['num_side_ss',num2str(i),'_dimID'];
        eval([name_elem_ss_varID '= netcdf.defVar(ncid,elem_ss_str,int_str,eval(num_side_ss_str));']);
        
        % set side_ss1, side_ss2, etc
        side_ss_str = ['side_ss',num2str(i)];
        name_side_ss_varID = genvarname(['side_ss',num2str(i),'_dimID']);
        eval([name_side_ss_varID '= netcdf.defVar(ncid,side_ss_str,int_str,eval(num_side_ss_str));']);
        
        % set dist_fact_ss1, dist_fact_ss2, etc
        dist_fact_ss_str = ['dist_fact_ss',num2str(i)];
        num_df_ss_str = ['num_df_ss',num2str(i)];
        name_dist_fact_ss_varID = genvarname(['dist_fact_ss',num2str(i),'_dimID']);
        eval([name_dist_fact_ss_varID '= netcdf.defVar(ncid,dist_fact_ss_str,double_str,eval(num_df_ss_str));']);
    end
end

elem_map_varID = netcdf.defVar(ncid,'elem_map','int',num_elem_dimID);
eb_status_varID = netcdf.defVar(ncid,'eb_status','int',num_el_blk_dimID);
eb_prop1_varID = netcdf.defVar(ncid,'eb_prop1','int',num_el_blk_dimID);
netcdf.putAtt(ncid,eb_prop1_varID,'name','ID')

for i =1:1:num_el_blk
    num_el_in_blk_str = ['num_el_in_blk',num2str(i),'_dimID'];
    num_nod_per_el_str = ['num_nod_per_el',num2str(i),'_dimID'];
    
    % set connect1, connect2, etc
    connect_str = ['connect',num2str(i)];
    name_connect_varID = genvarname(['connect',num2str(i),'_varID']);
    eval([name_connect_varID '= netcdf.defVar(ncid,connect_str,int_str,[eval(num_nod_per_el_str),eval(num_el_in_blk_str)]);']);
    switch eletype{i}
        case 'bar2'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','BAR2');
        case 'tri3'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','TRI3');
        case 'tri6'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','TRI6');
        case 'quad4'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','QUAD4');
        case 'tet4'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','TET4');
        case 'tet10'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','TET10');
        case 'quad8'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','QUAD8');
        case 'quad9'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','QUAD9');
        case 'hex8'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','HEX8');
        case 'hex20'
            netcdf.putAtt(ncid,eval(name_connect_varID),'elem_type','HEX20');
        otherwise
            error('incorrect element type');
    end
end

coord_varID = netcdf.defVar(ncid,'coord','double',[num_nodes_dimID,num_dim_dimID]);

% set global attibutes
% NOTE: some of the global attributes could not be set properly,
% instead they get copied from existing GENESIS file 'Default.g'
% Open Default.g file (should be under .../Matlab/misc/
ncid_tmp = netcdf.open([pwd '/misc/Default.g'],'NC_NOWRITE');
% copy values from Default.g
netcdf.copyAtt(ncid_tmp,netcdf.getConstant('NC_GLOBAL'),'api_version',ncid,netcdf.getConstant('NC_GLOBAL'));
netcdf.copyAtt(ncid_tmp,netcdf.getConstant('NC_GLOBAL'),'version',ncid,netcdf.getConstant('NC_GLOBAL'));
netcdf.copyAtt(ncid_tmp,netcdf.getConstant('NC_GLOBAL'),'floating_point_word_size',ncid,netcdf.getConstant('NC_GLOBAL'));
netcdf.copyAtt(ncid_tmp,netcdf.getConstant('NC_GLOBAL'),'file_size',ncid,netcdf.getConstant('NC_GLOBAL'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'file_size',1);

% close Default.g file
netcdf.close(ncid_tmp)

% set remaining global attributes
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title',title);

% leave 'Def'-mode
netcdf.endDef(ncid);

% write qa_records
netcdf.putVar(ncid,qa_records_varID,[0; 0; 0],[5; 1; 1],'CUBIT');
netcdf.putVar(ncid,qa_records_varID,[0; 1; 0],[4; 1; 1],'11.1');
netcdf.putVar(ncid,qa_records_varID,[0; 2; 0],[10; 1; 1],datestr(now, 'mm/dd/yyyy'));
netcdf.putVar(ncid,qa_records_varID,[0; 3; 0;],[8; 1; 1],datestr(now, 'HH:MM:SS'));

% write coor_names
for i = 1:1:length(coordinate_names)
    netcdf.putVar(ncid,coor_names_varID,[0; 0+i-1],coordinate_names{i});
end

for i = 1:1:num_node_sets
    % write ns_status
    netcdf.putVar(ncid,ns_status_varID,[0+i-1],ns_status(i));
    
    % write ss_prop1
    netcdf.putVar(ncid,ns_prop1_varID,[0+i-1],ns_prop1(i));
    
    % write node_ns1, node_ns2,etc
    node_ns_str = ['node_ns',num2str(i),'_dimID'];
    for j = 1:1:num_nod_ns(i) %%
        netcdf.putVar(ncid,eval(node_ns_str),[0+j-1],nodes_in_nodeset(i,j));
    end
    
    % write dist_fact_ns1, dist_fact_ns2,etc
    dist_fact_ns_str = ['dist_fact_ns',num2str(i),'_dimID'];
    for j = 1:1:num_df_ns(i)
        netcdf.putVar(ncid,eval(dist_fact_ns_str),[0+j-1],dist_fact_ns(i,j));
    end
    
end

for i = 1:1:num_side_sets
    % write ss_status
    netcdf.putVar(ncid,ss_status_varID,[0+i-1],ss_status(i));
    
    % write ss_prop1
    netcdf.putVar(ncid,ss_prop1_varID,[0+i-1],ss_prop1(i));
    
    % write elem_ss1, elem_ss2,etc
    elem_ss_str = ['elem_ss',num2str(i),'_dimID'];
    for j = 1:1:num_side_ss(i)
        netcdf.putVar(ncid,eval(elem_ss_str),[0+j-1],ele_in_sideset(i,j));
    end
    
    % write side_ss1, side_ss2,etc
    side_ss_str = ['side_ss',num2str(i),'_dimID'];
    for j = 1:1:num_side_ss(i)
        netcdf.putVar(ncid,eval(side_ss_str),[0+j-1],side_in_sideset(i,j));
    end
    
    % write dist_fact_ss1, dist_fact_ss2,etc
    dist_fact_ss_str = ['dist_fact_ss',num2str(i),'_dimID'];
    for j = 1:1:num_df_ss(i)
        netcdf.putVar(ncid,eval(dist_fact_ss_str),[0+j-1],dist_fact_ss(i,j));
    end
    
end

% write elem_map
for i = 1:1:numfele
    netcdf.putVar(ncid,elem_map_varID,[0+i-1],i);
end

for i = 1:1:num_el_blk
    % write eb_status
    netcdf.putVar(ncid,eb_status_varID,i-1,eb_status_val(i));
    % write eb_prop1
    eb_prop_str = ['eb_prop',num2str(i),'_dimID'];
    netcdf.putVar(ncid,eb_prop1_varID,i-1,eb_prop1_val(i));
end

% write coord
netcdf.putVar(ncid,coord_varID,[0 0]',[numfenod num_dim]',fecoord);
% for i = 1:1:numfenod
%     for d = 1:1:num_dim
%         netcdf.putVar(ncid,coord_varID,[i-1 d-1]',[1 1]',fecoord(i,d));
%     end
% end

% write connect1, connect2, etc.
ele_count = 1;
% ele_count = 0;

for i = 1:1:num_el_blk
    
    connect_str = ['connect',num2str(i),'_varID'];
    
    netcdf.putVar(ncid,eval(connect_str),[0 0]',[num_nod_per_el(i) num_el_in_blk(i)]',fetopo(ele_count:ele_count+num_el_in_blk(i)-1,2:num_nod_per_el(i)+1)');
    ele_count = ele_count+num_el_in_blk(i);
    
%     for j = 1:1:num_el_in_blk(i)
%         ele_count = ele_count+1;
% 
%         for n = 1:1:num_nod_per_el(i)
%             netcdf.putVar(ncid,eval(connect_str),[n-1 j-1]',[1 1]',fetopo(ele_count,n+1));
%         end
%         
%     end
end

% close exodus file
netcdf.close(ncid)
