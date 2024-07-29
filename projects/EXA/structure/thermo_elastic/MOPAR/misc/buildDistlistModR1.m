%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% REMEMBER TO DELETE LAST EMPTY LINE FROM DISTTABLE FILE!!!
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input

% Define variables, folderName without / at the end
 folderName='/home/tkachuk/codes/Level_set/phi_2D_for_3D_half_conv_concave_sym_filter';
%folderName='/home/tkachuk/Apps/cubit.11.1';
meshName='phi_2D_for_3D.g';
rmax=1.6; %0.01;  % Smoothing radius
numlev=4;   % Levels of adjacent elements to search for nodes in smoothing radius 
nproc=1;    % Number of processors to use in the search

saveFile = 1;

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

xplane = 0; 
yplane = NaN;

tol = 1e-9;
% build ghost nodes and "full" node list 

cNumNodes = numfenod;
% put original nodes into new list
node(numfenod).coords = [0,0,0];
node(numfenod).elems = [];
node(numfenod).origInd = [0];
for i = 1:numfenod
    node(i).coords = fecoord(i,:);
    node(i).elems = find(sum(ismember(fetopo(:,2:nnpe),i),2)>0);
    node(i).origInd = i;
end

% try to do x plane reflection
if ~isnan(xplane)
        for i = 1:numfenod
            distToPlane = norm(node(i).coords(1) - xplane);
            
            if distToPlane < tol
                distToPlane = inf;
            end
        
            if distToPlane <= rmax
                newX = 2*xplane - node(i).coords(1);
                
                cNumNodes = cNumNodes + 1;
                node(cNumNodes).coords = [newX,node(i).coords(2:3)];
                node(cNumNodes).elems = node(i).elems;
                node(cNumNodes).origInd = node(i).origInd;
            end
        end
end

numNodesBeforeYFlip = cNumNodes;
if ~isnan(yplane)
    for i = 1:numNodesBeforeYFlip
        distToPlane = norm(node(i).coords(2) - yplane);
        
        if distToPlane < tol
            distToPlane = inf;
        end
        
        if distToPlane <= rmax
            newY = 2*yplane - node(i).coords(2);

            cNumNodes = cNumNodes + 1;
            node(cNumNodes).coords = [node(i).coords(1),newY,node(i).coords(3)];
            node(cNumNodes).elems = node(i).elems;
            
            node(cNumNodes).origInd = node(i).origInd;
        end
    end
end

numTotalNodes = cNumNodes;

maxNumElem = length(fetopo(:,1));

for i = 1:maxNumElem
    elem(i).topo = fetopo(i,2:end);
end

for i = (numfenod+1):numTotalNodes
    myElems = node(i).elems;
    for j = 1:length(myElems)
        elem(myElems(j)).topo = [elem(myElems(j)).topo,i];
    end
end

% elem(maxNumElem).topo = [];
% elemNodes = zeros(maxNumElem,1);

% maxNodePerElem = 0;
% for i = 1:numTotalNodes
%     for j = 1:length(node(i).elems)
%         
%         if (elemNodes(node(i).elems(j)) == 0) || (~any(elem(node(i).elems(j)).topo(1:elemNodes(node(i).elems(j))) == i))
%             elemNodes(node(i).elems(j)) = elemNodes(node(i).elems(j)) + 1;
%             elem(node(i).elems(j)).topo(elemNodes(node(i).elems(j))) =  i;
%             
%             node(i).elems(j)
%             if any(elem(node(i).elems(j)).topo > 60)
%                 disp('ack');
%             end
%             
%             if length(elem(node(i).elems(j)).topo) > maxNodePerElem
%                 maxNodePerElem = length(elem(node(i).elems(j)).topo) ;
%             end
%         end
%     end
% end

fetopoNew = zeros(maxNumElem,maxNodePerElem);
for i = 1:maxNumElem
    nNodesThisElem = length(elem(i).topo);
    fetopoNew(i,1:nNodesThisElem) = elem(i).topo;
end

fetopoNew = [(1:maxNumElem)',fetopoNew];

fecoordNew = zeros(numTotalNodes,3);
nodeInd = zeros(numTotalNodes,1);
% figure
for i  = 1:numTotalNodes
    fecoordNew(i,:) = node(i).coords;
    nodeInd(i) = node(i).origInd;
%     plot3(node(i).coords(1),node(i).coords(2),node(i).coords(3),'.');
%     hold on
end
% grid on

% pause

% compute distances

parfor ip=1:nproc
    
    fname=sprintf([folderName,'table_%03d'],ip);

    if saveFile
        fid=fopen(fname,'w');
    end
        
    lfecoord=fecoordNew;
    lfetopo =fetopoNew;
    
    for ik=1:nnpp(ip)
        
        in=lnps(ip,ik);
        
        nodlist=in;
        
        for ilev=1:numlev
            [faa]=sum(ismember(lfetopo(:,2:end),nodlist),2)>0;
            nodlist=unique(lfetopo(faa,2:end));
            nodlist = sort(nodlist);
            if nodlist(1) == 0
                nodlist = nodlist(2:end);
            end
        end
        
        nndl = length(nodlist);
        ww   = zeros(nndl,1);
        
        for ii=1:nndl
            dist=norm(lfecoord(nodlist(ii),:)-fecoordNew(in,:));
            dels=rmax-dist;
            if dels < 0; ww(ii)=-1; else ww(ii)=dels; end
        end
        
        nodlist(ww==-1)=[];
        ww(ww==-1)=[];
        nndl = length(nodlist);
        
         ww=1/sum(ww)*ww; 
         
         % combine ghost with their original nodes
         nodlistNew = unique(nodeInd(nodlist));
         nndlNew = length(nodlistNew);
         wwNew = zeros(1,nndlNew);
         
         for myI = 1:nndl
             for myJ = 1:nndlNew
                 if nodlistNew(myJ) == nodeInd(nodlist(myI))
                     wwNew(myJ) = wwNew(myJ) + ww(myI);
                 end
             end
         end
         
         if saveFile
            fprintf(fid,'%1.7e %1.7e ',in,nndlNew);
            for ii=1:nndlNew
                fprintf(fid,'%1.7e ',nodlistNew(ii));
            end
            for ii=1:nndlNew
                fprintf(fid,'%1.7e ',wwNew(ii));
            end

            if ~((ip==nproc)&&(ik==nnpp(ip)))
                fprintf(fid,'\n');
            end
         end

    end
    
    if saveFile
        fclose(fid);
    end
    
end

if saveFile
    cmd=['rm -f ',folderName,'disttable'];
    system(cmd);
    for ip=1:nproc
        cmd=sprintf(['cat ',folderName,'table_%03d >> ',folderName,'disttable'],ip);
        system(cmd);
    end
end
