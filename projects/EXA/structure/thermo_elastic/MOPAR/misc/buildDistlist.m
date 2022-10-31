%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% REMEMBER TO DELETE LAST EMPTY LINE FROM DISTTABLE FILE!!!
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input

% Define variables, folderName without / at the end
folderName='/home/maute';
meshName='FullPlate.g';
rmax=0.6;  % Smoothing radius
numlev=4;   % Levels of adjacent elements to search for nodes in smoothing radius 
nproc=6;    % Number of processors to use in the search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add backslash to folder name
folderName=[folderName,'/'];

% Append folder name to mesh name
meshName=[folderName,meshName];

% Delete previous tables
cmd=['rm -f ',folderName,'table_*'];
system(cmd);

% Read mesh
[fecoord,fetopo,~,~,numfenod,numfele] = LoadExoMesh(meshName);

% Decompose with respect to nodes
nnpe = size(fetopo,2);
nnps=ceil(numfenod/nproc);
nnpp=zeros(nproc,1);
lnps=zeros(nproc,nnps);
ins=0;

for j=1:nnps
    for i=1:nproc
        ins=ins+1; if ins <= numfenod;nnpp(i)=nnpp(i)+1; end
    end
end

if sum(nnpp) ~= numfenod;error('incorrect node-processor assignment'); end

ins=0; ies=0; ics=0; ibs=0;

for i=1:nproc
    inx=ins+nnpp(i);  
    lnps(i,1:nnpp(i))=ins+1:inx;
    ins=inx;
end

% compute distances

parfor ip=1:nproc
    
    fname=sprintf([folderName,'table_%03d'],ip);

    fid=fopen(fname,'w');
        
    lfecoord=fecoord;
    lfetopo =fetopo;
    
    for ik=1:nnpp(ip)
        
        in=lnps(ip,ik);
        
        nodlist=in;
        
        for ilev=1:numlev
            [faa]=sum(ismember(lfetopo(:,2:nnpe),nodlist),2)>0;
            nodlist=unique(lfetopo(faa,2:nnpe));
        end
        
        nndl = length(nodlist);
        ww   = zeros(nndl,1);
        
        for ii=1:nndl
            dist=norm(lfecoord(nodlist(ii),:)-fecoord(in,:));
            dels=rmax-dist;
            if dels < 0; ww(ii)=-1; else ww(ii)=dels; end
        end
        
        nodlist(ww==-1)=[];
        ww(ww==-1)=[];
        nndl = length(nodlist);
        
         ww=1/sum(ww)*ww; 
         
        fprintf(fid,'%1.7e %1.7e ',in,nndl);
        for ii=1:nndl
            fprintf(fid,'%1.7e ',nodlist(ii));
        end
        for ii=1:nndl
            fprintf(fid,'%1.7e ',ww(ii));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
    
end

cmd=['rm -f ',folderName,'disttable'];
system(cmd);
for ip=1:nproc
    cmd=sprintf(['cat ',folderName,'table_%03d >> ',folderName,'disttable'],ip);
    system(cmd);
end
