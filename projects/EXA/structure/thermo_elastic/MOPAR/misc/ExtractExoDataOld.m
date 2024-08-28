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
    elseif strcmp(text(1:12),'name_nod_var');
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
    elseif strcmp(text(1:13),'name_elem_var');
        EVnames = netcdf.getVar(ncid,iv);
        varstring = sprintf('EVar{}.Name');
    elseif strcmp(text(1:8),'elem_map');
        elem_map=netcdf.getVar(ncid,iv);
        varstring = sprintf('GVar{}.Name');
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
if ~isempty(Coordnames);
    for ic = 1:ci
        Coord{ic}.Name = Coordnames(:,ic)';
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

% if specified save the data in a matlab file
if nargin == 3
    save(savefile,'Data')
end
