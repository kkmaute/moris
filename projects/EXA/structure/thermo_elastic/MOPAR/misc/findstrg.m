function line=awk(filename,pattern,columid,delim,numflag)

if nargin<4
  delim=' ';
end

if nargin<5
    numflag=0;
end

fid = fopen(filename,'r');
C = textscan(fid, '%s','Delimiter','\n');
fclose(fid);
C = C{:};

Lia = ~cellfun(@isempty, strfind(C,pattern));

line = cell(sum(Lia),1);

k=0;
for i=1:length(Lia)
    if Lia(i)
        k=k+1;
        if nargin>2
            abc=strsplit(C{i},delim);
            line{k}=abc{columid};
        else
            line{k} = C{i};
        end
    end
end

if numflag
    nvalue=zeros(length(line),1);
    for i=1:length(line)
        nvalue(i)=str2double(line{i});
    end
    line=nvalue;
end
