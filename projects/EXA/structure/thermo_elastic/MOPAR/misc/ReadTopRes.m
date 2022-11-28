function  [res,times,key]=ReadTopRes(resfile)

maxtimes=100;

fid=fopen(resfile,'r');

% read first line
line=fgetl(fid);

% read output type and identifier
ids=textscan(line,'%s');

stp = ids{1}(1);
key = ids{1}(2);

type=0;
if strcmp('Scalar',stp); type=1;end
if strcmp('Vector',stp); type=2;end

switch type
    
    case 1
        resdim=1;
    case 2
        resdim=3;
    otherwise
        error('incorrect data type');
end

% read number of nodes
line=fgetl(fid);
numnodes=sscanf(line,'%d');

it=0;

times=zeros(maxtimes,1);
res=cell(maxtimes,1);

while 1
    
    line=fgetl(fid);
    
    if line < 0;break;end
   
    it=it+1;
    
    times(it)=sscanf(line,'%e');
    
    tres=zeros(numnodes,resdim);
    
    for in=1:numnodes
        
        line=fgetl(fid);
        
        if line < 0
            fprintf(' ... incomplete time set\n');
            it=it-1;
            break;
        end
        
        switch type
            
            case 1
                tres(in,:)=sscanf(line,'%e');
            case 2
                tres(in,:)=sscanf(line,'%e %e %e');
            otherwise
                error('incorrect data type');
        end
    end

    res{it}=tres;
    
end

fclose(fid);

times=times(1:it);
res=res(1:it);

