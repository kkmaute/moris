function ResetTimeExoDataForESsuffix(exofile)
%
% adjusts global time variable in sequence of exofiles
%

cmd=sprintf('ls %s* | wc',exofile);

[irc,wcstr]=system(cmd);
wcnum=str2num(wcstr);
wcnum=323;

Time=0;

for i=172:wcnum(1)
    
    xfemxfn=sprintf('device.e-s.%04d',i); %[exofile '.e-s.000' num2str(i)];
    
    if ~exist(xfemxfn); continue; end
    
    % open exodus file
    ncid = netcdf.open(xfemxfn,'NC_WRITE');
    
    Time=Time+1;      % global time
    TimeStep=1;       % time step in file 
                      % (assuming that there is only one time step per file)
    
    varid = netcdf.inqVarID(ncid,'time_whole');
    oldtime = netcdf.getVar(ncid,varid,'double');
    netcdf.putVar(ncid,varid,TimeStep-1,Time);
    newtime = netcdf.getVar(ncid,varid,'double');
    
    fprintf('%s   %f -> %f\n',xfemxfn,oldtime,newtime);
    
    % close exodus file
    netcdf.close(ncid)
    
end
