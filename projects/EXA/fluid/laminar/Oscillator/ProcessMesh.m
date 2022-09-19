%
% script to write NASTRAN input file as exodus mesh 
% considering multiple blocks but NO side set information
%

clc;
close all;
clear all;

% name of output exodus file
outputfile='Oscillator.exo';


% tools for reading and writing exodus meshes
scripts_dir = './misc';
addpath(scripts_dir);

% read fe data
%[Fedata,Eledata_tri6,Eledata_quad8]=elefenasdata();
%[Fedata]=elefenasdata_1();
[Fedata, Eledata_tri6,Eledata_quad8] = DataExtractor();
%[Eledata_tri6,Eledata_quad8]=Input_Data();

% echeck that node numbers are consecutive
if max(Fedata(:,1))~=size(Fedata,1)
    error('node numbers not consecutive - feature not implemented')
end

% process data
numfenod=size(Fedata,1);
fecoord=Fedata(:,2:3);

numtri6=size(Eledata_tri6,1);
numquad8=size(Eledata_quad8,1);

numfele=numtri6+numquad8;

fetopo=zeros(numfele,11);

fetopo(1:numtri6,1:2)=Eledata_tri6(:,1:2);
fetopo(1:numtri6,3)=6;
fetopo(1:numtri6,4:9)=Eledata_tri6(:,3:end);

fetopo(numtri6+1:end,1:2)=Eledata_quad8(:,1:2);
fetopo(numtri6+1:end,3)=8;
fetopo(numtri6+1:end,4:11)=Eledata_quad8(:,3:end);

numblks_tri6=max(Eledata_tri6(:,2));
numblks_quad8=max(Eledata_quad8(:,2));
numblks=numblks_tri6+numblks_quad8;

fetopo(numtri6+1:end,2)=fetopo(numtri6+1:end,2)+numblks_tri6;

eletypeidx=zeros(numblks,1);

for ie=1:numfele
    blk=fetopo(ie,2);
    typdix=fetopo(ie,3);
    
    if typdix==6;typdix=3;end
    if typdix==8;typdix=4;end
    
    if eletypeidx(blk) == 0
        eletypeidx(blk)=typdix;
    else
        if eletypeidx(blk) ~= typdix
            error('inconsistent element index');
        end
    end
end

eletypes=cell(numblks,1);

for ib=1:numblks
    switch eletypeidx(ib)
        case 3
            eletypes{ib}='tri3';
        case 6
            eletypes{ib}='tri6';
        case 4
            eletypes{ib}='quad4';
        case 8
            eletypes{ib}='quad8';
        otherwise
            error('element type not recognized');
    end
end

cordname={'x','y'};

% just consider volume elements
numvolele=numfele; %sum(fetopo(:,1)==0);

numeleinblks=zeros(numblks,1);
for ie=1:numfele
    blk=fetopo(ie,2);
    numeleinblks(blk)=numeleinblks(blk)+1;
end

fetopoout=zeros(numvolele,size(fetopo,2)-2);
eleinblkscounter=zeros(numblks,1);

for ie=1:numfele
    blk=fetopo(ie,2);
    eleidx=eletypeidx(blk);
    eleinblkscounter(blk)=eleinblkscounter(blk)+1;
%     fedata(ie)=femdata.EVar{6,blk}.Val(eleinblkscounter(blk));
    neweleidx=eleinblkscounter(blk)+sum(numeleinblks(1:blk-1));
    fetopoout(neweleidx,1:eleidx+1)=[eleidx fetopo(ie,4:3+eleidx)];
end

CreateNetcdfFile(outputfile,eletypes,cordname,fecoord,fetopoout, ...
    numvolele,zeros(numfenod,1),numblks,numeleinblks);
IniExoData(outputfile,fecoord);

% comment out to visualize side sets
%
% sidemaptri=[1,2;2,3;3,1];
% sidemapquad=[1,2;2,3;3,4;4,1];
% 
% for ie=1:numfele
%     if fetopo(ie,1)~=0
%         numnodeperele=fetopo(ie,3);
%         sideordinal=feprop(ie,1);
%         switch numnodeperele
%             case 3
%                 nodeidx=fetopo(ie,4:6);
%                 sidenodes=nodeidx(sidemaptri(sideordinal,:));
%             case 4
%                 nodeidx=fetopo(ie,4:7);
%                 sidenodes=nodeidx(sidemapquad(sideordinal,:));
%             otherwise
%                 error('strange')
%         end
%         plot(fecoord(sidenodes,1),fecoord(sidenodes,2),'r.'); hold on
%     end
% end
