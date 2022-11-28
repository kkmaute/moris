function ResetTimeExoData(basename,tstepsize,wcnum,sourcefile)
%
% adjusts global time variable in sequence of exofiles
%

Time=-tstepsize;

copyflag=0;
if nargin>3
    copyflag=1;
end

i=0;

while i<wcnum
    
    i=i+1;
    
    xfemxfn=sprintf('%s.e-s.%04d',basename,i); %[exofile '.e-s.000' num2str(i)];
    
    if copyflag && i>1
        copyfile(sourcefile,xfemxfn);
    end

    fprintf('processing %s\n',xfemxfn);
         
    if ~exist(xfemxfn); continue; end
    
    % open exodus file
    ncid = netcdf.open(xfemxfn,'NC_WRITE');
    
    Time=Time+tstepsize;      % global time
    TimeStep=1;               % time step in file 
                              % (assuming that there is only one time step per file)
    
    varid = netcdf.inqVarID(ncid,'time_whole');
    oldtime = netcdf.getVar(ncid,varid,'double');
    netcdf.putVar(ncid,varid,TimeStep-1,Time);
    newtime = netcdf.getVar(ncid,varid,'double');
    
    fprintf('%s   %f -> %f\n',xfemxfn,oldtime,newtime);
    
    % close exodus file
    netcdf.close(ncid)
    
end
