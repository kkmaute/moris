function top2exo(topfile,resfile,exofile)
%
% conversion of mesh and results in top format to exodus
%
% input : topfile - name of file with mesh in top formation (string)
%         resfile - cell array with name of result files in top formation (cell arrray of strings)
%         exofile - name of exodus output file (string)
%

fprintf(' ... read mesh file\n');

[fecoord,fetopo,neib,ebtyp]=ReadTopMesh(topfile);

fprintf(' ... create exo mesh\n');

CreateNetcdfFile(exofile,ebtyp,{'x','y','z'},fecoord,fetopo,size(fetopo,1),ones(size(fecoord,1),1),size(neib,1),neib);

IniExoData(exofile,fecoord);

fprintf(' ... read result file(s)\n');

sufix=['x' 'y' 'z'];
ik=0;

for ir=1:length(resfile)
    
    [res,times,key]=ReadTopRes(resfile{ir});
    
    resdim=size(res{1},2);
    
    if resdim == 1
        NodVarName{ik+1}=key{1};
    else
        for id=1:resdim
            NodVarName{ik+id}=[key{1} '' sufix(id)];
        end
    end
    
    for it=1:length(times)
        for id=1:resdim
            resvec=res{it}(:,id);
            allres{it,ik+id}=resvec;
        end
    end
    ik=ik+resdim;
end

fprintf(' ... create exo results\n');

for it=1:length(times)
    WriteExoData(exofile,NodVarName,[allres{it,:}],it,times(it));
end
