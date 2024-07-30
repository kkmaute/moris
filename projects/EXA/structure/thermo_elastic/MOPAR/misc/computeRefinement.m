
clear all
close all
clc

datafile ='bend.g'; %'bend.e-s.0300'; %'bend_coarse.exo';
levsbound= 0.003;
refstrat = 3;
% read data

Data = ExtractExoData(0,datafile);

% extract information

fecoord=[ Data.Coord{1}.Coords Data.Coord{2}.Coords Data.Coord{3}.Coords ];

topo=Data.Etopo{1}.Nodes;

nlevs=Data.NVar{5}.Val(:,1);

numfenod=size(nlevs,1);
numfele=size(topo,2);

% compute elemental nodal level set values for blockset 1

eleFlag  = 0;

switch refstrat
    case 1 % based on elemental levelset value

        eleFlag = 1;
        
        elevs = zeros(numfele,1);
        
        for ie=1:numfele
            elevs(ie)=mean(nlevs(topo(:,ie)));
        end

        tobeRefined=intersect(find(elevs > -levsbound),find(elevs < levsbound));

    case 2 % based on noal levelset value        
        
        usedn = zeros(numfenod,1);

        for ie=1:numfele
            usedn(topo(:,ie))=1;
        end
        
        nlevs(usedn==0)=10*levsbound;
        
        tobeRefined=intersect(find(nlevs > -levsbound),find(nlevs < levsbound));
    
    case 3 % based on whether elements are intersected
        
        eleFlag = 1;

        isect = zeros(numfele,1);

        for ie=1:numfele
            elevs=nlevs(topo(:,ie));
            isect(ie)=min(elevs)*max(elevs)<0;
        end
        
        tobeRefined=find(isect);
    
    case 4 % nodes belonging to intersected elements
        
        isectn = zeros(numfenod,1);

        for ie=1:numfele
            elevs=nlevs(topo(:,ie));
            if min(elevs)*max(elevs)<0
                isectn(topo(:,ie))=1;
            end
        end
        
        tobeRefined=find(isectn);
        
    otherwise
        error('incorrect refinement strategy');
end        
        
fid=fopen('refineMesh.txt','w');

if (eleFlag)
    fprintf('number of elements to be refined: %d\n\n',length(tobeRefined));
    fprintf(fid,'set duplicate block elements on\nblock 111 hex ');
else
    fprintf('number of nodes to be refined: %d\n\n',length(tobeRefined));
    fprintf(fid,'nodeset 111 node ');
end

for ir=1:length(tobeRefined)
    fprintf(fid,'%d,',tobeRefined(ir));
end

if (eleFlag)
    fprintf(fid,'\nrefine hex ');
else
    fprintf(fid,'\nrefine node ');
end

for ir=1:length(tobeRefined)
    fprintf(fid,'%d,',tobeRefined(ir));
end

fprintf(fid,' numsplit 1 bias 1 depth 1 smooth');

fclose(fid);

