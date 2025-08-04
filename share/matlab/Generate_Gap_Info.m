function Generate_Gap_Info

clear all
close all

global nlbgsiter newiter lastemptydata lclast

newiter=0;          % newton iteration count for nlbgs
nlbgsiter=1;        % nlbgs iteration count for nlbgs
lastemptydata=100;  % last time empty data was written
lclast=1;           % last line number that as processed

%read logfile
logfile = readlines('log');

keywordlist={'Mlika','Converged','NotConverged','NormalMlika','TractionMlika'}';

nlbgs=1;          % flag for nlbgs: 1 nlbgs used, 0 only newton used
keyid=1;          % keyword id
maxiter=669;      % maximum number of iterations
delp=-50;          % range of pressure values (for scaling)
delaug=-50;        % range of aug term values (for scaling)
deltract=-50;      % range of traction magnitude (for scaling)

for iter=1:maxiter
    stopflag=read_newton_iteration(logfile,nlbgs,iter-1,keywordlist,keyid,delp,delaug,deltract);
    if stopflag
        break;
    end
end
end

%============================================
function stopflag=read_newton_iteration(logfile,nlbgs,iter,keywordlist,keyid,delp,delaug,deltract)

global lastemptydata

keyword=keywordlist{keyid};

data=process_log_file(logfile,nlbgs,iter,keyword);

stopflag=0;
if isempty(data)
    if lastemptydata<iter-1
        stopflag=1;
        return;
    end
    lastemptydata=iter;
end

if keyid<4
    read_gap(iter,data,keyword,delp,delaug);
else
    read_normal_traction(iter,data,keyword,deltract);
end

end

%============================================

function data=process_log_file(logfile,nlbgs,iter,keyword)

global nlbgsiter newiter lclast

tmpdata={};
data=[];

lcmax=size(logfile,1);

if nlbgs 
    % read until we find: NLBGS - Iteration: nlbgsiter
    lc=lclast;
    keystr=sprintf('NLBGS - Iteration: %d',nlbgsiter);
    while lc<lcmax && contains(logfile(lc),keystr) == 0 
        lc=lc+1;
    end
    if lc==lcmax
        return;
    end
    
    % read until we find: NLBGS - Iteration: nlbgsiter + 1 and find maximum
    % Newton iteration
    kc=lc+1;
    maxnewiter=-1;
    keystr=sprintf('NLBGS - Iteration: %d',nlbgsiter+1);
    while  kc<lcmax && contains(logfile(kc),keystr) == 0
        if contains(logfile(kc),'Newton - Iteration:') == 1
            parts=split(logfile(kc),'Newton - Iteration:');
            maxnewiter=max(maxnewiter,str2double(parts(end)));
        end
        kc=kc+1;
    end
    if maxnewiter<0
        return;
    end
    
    % update iteration counters
    newiter=newiter+1;
    if newiter>maxnewiter
        nlbgsiter=nlbgsiter+1;
        newiter=1;
        lc=kc; 
        kc=lc+1;
        keystr=sprintf('NLBGS - Iteration: %d',nlbgsiter+1);
        while kc<lcmax && contains(logfile(kc),keystr) == 0
            kc=kc+1;
        end
        kc=min(kc,lcmax);
    end
    lclast=lc;

    % scan log file between lines lc and kc for Niter = newiter keyword"
    keystr=sprintf('Niter = %d %s',newiter,keyword);
    for i=lc:kc
        if contains(logfile(i),keystr) == 1
            parts=split(logfile(i),keyword);
            tmpdata{end+1}=str2num(parts(end));
        end
    end
else
    keystr=sprintf('Niter = %d %s',iter,keyword);
    keystrnext=sprintf('Niter = %d %s',iter+1,keyword);
    for i=lclast:lcmax
        if contains(logfile(i),keystr) == 1
            parts=split(logfile(i),keyword);
            tmpdata{end+1}=str2num(parts(end));
        end
        % check if we reached next iteration
        % if yes, store line number and exit
        if contains(logfile(i),keystrnext) == 1
            lclast=i;
            break;
        end
    end
end

%convert to matrix
if ~isempty(tmpdata)
    numcols=size(tmpdata{1},2);
    data=zeros(length(tmpdata),numcols);
    for i=1:size(data,1)
        data(i,:)=tmpdata{i};
    end
end
end

%============================================

function read_gap(iter,data,keyword,delp,delaug)

% number of normals
numNormals=size(data,1);

fprintf('Iteration: %d number of normals: %d\n',iter+1,numNormals);

if numNormals == 0
    data=[0 0 0 0];
    numNormals=1;
end

% number of points
numPoints=numNormals*2;

delx1=max(data(:,1))-min(data(:,1));
dely1=max(data(:,2))-min(data(:,2));
delx2=max(data(:,3))-min(data(:,3));
dely2=max(data(:,4))-min(data(:,4));
maxx=max([delx1,delx2,dely1,dely2]);

pscale=0.0;
augscale=0.0;

if size(data,2) > 4
    pressure=data(:,5);
    augterm=data(:,6);
    gap=data(:,7);
    pmax=max(data(:,5));
    pmin=min(data(:,5));
    if delp<0
        delp=pmax-pmin;
    end
    if delp>1e-9
        pscale=0.1*maxx/delp;
    end
    augmax=max(data(:,6));
    augmin=min(data(:,6));
    if delaug<0
        delaug=augmax-augmin;
    end
    if delaug>1e-9
        augscale=0.1*maxx/delaug;
    end
    gmax=max(data(:,7));
    gmin=min(data(:,7));

    fprintf('pressure range: min = %e max = %e  delta = %e  pscale   = %e\n',pmin,pmax,delp,pscale);
    fprintf('Aug term range: min = %e max = %e  delta = %e  augscale = %e\n',augmin,augmax,delaug,augscale);
    fprintf('gap range:      min = %e max = %e  delta = %e  \n',gmin,gmax,gmax-gmin);
else
    pressure=zeros(numNormals,1);
    augterm=zeros(numNormals,1);
    gap=zeros(numNormals,1);
end

points=zeros(numPoints,3);
normals=zeros(numNormals,2);
inod=1;

for in=1:numNormals
    points(inod,1:2)=data(in,1:2);
    points(inod,3)=pscale*pressure(in);
    inod=inod+1;
    points(inod,1:2)=data(in,3:4);
    points(inod,3)=pscale*pressure(in);
    inod=inod+1;
    normals(in,:)=[inod-3 inod-2];  %uses null based index
end

% Open file
filename=sprintf('normals_%s_%03d.vtk',keyword,iter+1);
cmd=sprintf('rm -f %s',filename);
system(cmd);

fid = fopen(filename, 'w');

% Write header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from MATLAB\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');

% Write points
fprintf(fid, 'POINTS %d float\n', size(points, 1));
fprintf(fid, '%e %e %e\n', points');

% Write vetices
fprintf(fid, 'VERTICES %d %d\n',numNormals, 2*numNormals);
for i = 1:numNormals
    fprintf(fid, '1 %d \n', normals(i, 1));
end

% Write lines
fprintf(fid, 'LINES %d %d\n', numNormals, 3*numNormals );
for i = 1:numNormals
    fprintf(fid, '2 %d %d\n', normals(i, 1), normals(i, 2));
end

% Write scalar data for lines
fprintf(fid, 'CELL_DATA %d\n', 2*size(pressure, 1));
fprintf(fid, 'SCALARS pressure float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%e\n', pressure);
fprintf(fid, '%e\n', pressure);

fprintf(fid, 'SCALARS augterm float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%e\n', augterm);
fprintf(fid, '%e\n', augterm);

fprintf(fid, 'SCALARS gap float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%e\n', gap);
fprintf(fid, '%e\n', gap);

% Close file
fclose(fid);
end

%============================================

function read_normal_traction(iter,data,keyword,deltract)

% number of normals
numNormals=size(data,1);

fprintf('Iteration: %d number of normals: %d\n',iter+1,numNormals);

if numNormals == 0
    data=[0 0 0 0 0 0];
    numNormals=1;
end

delx=max(data(:,1))-min(data(:,1));
dely=max(data(:,2))-min(data(:,2));
maxx=max([delx,dely]);

vec1scale=0.0;
vec2scale=0.0;

maxvec1=-1;
maxvec2=-1;
for in=1:numNormals
    maxvec1=max(maxvec1,norm(data(in,3:4)));
    maxvec2=max(maxvec2,norm(data(in,5:6)));
end

if deltract<0
    vec1scale=maxvec1;
    vec2scale=maxvec2;
end
if maxvec1>1e-9
    vec1scale=0.1*maxx/maxvec1;
end
if maxvec2>1e-9
    vec2scale=0.1*maxx/maxvec2;
end
fprintf('vector 1 range: %f  vec1scale: %f\n',maxvec1,vec1scale);
fprintf('vector 2 range: %f  vec2scale: %f\n',maxvec2,vec2scale);

% number of points
numPoints=numNormals*2;

points=zeros(numPoints,3);
normals=zeros(2*numNormals,2);

inod=1;
for in=1:numNormals
    points(inod,1:2)=data(in,1:2);
    inod=inod+1;
    points(inod,1:2)=data(in,1:2)+vec1scale*data(in,3:4);
    inod=inod+1;
    normals(in,:)=[inod-3 inod-2];  %uses null based index
end
for in=1:numNormals
    points(inod,1:2)=data(in,1:2);
    inod=inod+1;
    points(inod,1:2)=data(in,1:2)+vec2scale*data(in,5:6);
    inod=inod+1;
    normals(in+numNormals,:)=[inod-3 inod-2];  %uses null based index
end

% Open file
filename=sprintf('normals_%s_%03d.vtk',keyword,iter+1);
cmd=sprintf('rm -f %s',filename);
system(cmd);

fid = fopen(filename, 'w');

% Write header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from MATLAB\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');

% Write points
fprintf(fid, 'POINTS %d float\n', size(points, 1));
fprintf(fid, '%e %e %e\n', points');

% Write lines
fprintf(fid, 'LINES %d %d\n', 2*numNormals, 6*numNormals );
for i = 1:2*numNormals
    fprintf(fid, '2 %d %d\n', normals(i, 1), normals(i, 2));
end

% Write scalar data for lines
fprintf(fid, 'CELL_DATA %d\n', 2*numNormals);
fprintf(fid, 'SCALARS vectype float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%e\n', ones(numNormals,1));
fprintf(fid, '%e\n', 2*ones(numNormals,1));

% Close file
fclose(fid);
end
