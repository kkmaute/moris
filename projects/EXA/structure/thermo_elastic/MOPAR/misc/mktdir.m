function dname=mktdir(name,tpath)

if nargin == 0
    dname='.';
    return
end

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

dname=sprintf('%s/%s_%.5d',tpath,name,round(10000*rand));

mkdir(dname);
