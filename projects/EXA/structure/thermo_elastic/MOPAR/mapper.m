inputmesh='allmat.e';
outputmesh='exomesh.e';

exofile='test.e';

%==========================================================================
% set path to source directory
scripts_dir = './misc';
addpath(scripts_dir);


[inpfecoord,inpfetopo,inpfeprop,inpfebc,inpnumfenod,inpnumfele]= LoadExoMesh(inputmesh,1);
[outfecoord,outfetopo,outfeprop,outfebc,outnumfenod,outnumfele]= LoadExoMesh(outputmesh,2);

InpData=ExtractExoData(0,inputmesh);

inpTemp=InpData.NVar{1}.Val;

if inpnumfenod~=outnumfenod
    error('node numbers do not match');
end

innoddim=size(inpfecoord,2);
outnoddim=size(outfecoord,2);

noddim=min(innoddim,outnoddim);

nodiomap=zeros(inpnumfenod,1);

for in=1:inpnumfenod
    
    found=0;
    for io=1:outnumfenod
        
       if norm(inpfecoord(in,1:noddim)-outfecoord(io,1:noddim)) < 1e-6
           nodiomap(in)=io;
           found=1;
           break;
       end
    end
    if found==0
        error('node not found');
    end
end

% map temperatures
outTemp(nodiomap)=inpTemp;

cordname={'x','y'};
CreateNetcdfFile(exofile,{'tri3'},cordname,outfecoord,outfetopo, ...
    outnumfele,zeros(inpnumfenod,1),1,outnumfele);
IniExoData(exofile,outfecoord);

WriteExoData(exofile,{'Temp'},outTemp',1.0,0);

WriteExoData(outputmesh,{'Temp'},outTemp',1.0,0);


 