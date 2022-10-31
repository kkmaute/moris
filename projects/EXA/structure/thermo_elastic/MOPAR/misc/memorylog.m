% monitor memory use

printlog(sprintf('\n\n Global workspace (xfemfunc)\n\n'));
wsp=whos('global');
gbytes=0;
for i=1:length(wsp)
    gbytes=gbytes+wsp(i).bytes;
    printlog(sprintf('%s  = %e KB\n',wsp(i).name,wsp(i).bytes/1e3));
end

printlog(sprintf('\n\n Local workspace (xfemfunc)\n\n'));
wsp=whos;
tbytes=0;
for i=1:length(wsp)
    tbytes=tbytes+wsp(i).bytes;
    printlog(sprintf('%s  = %e KB\n',wsp(i).name,wsp(i).bytes/1e3));
end

printlog(sprintf('\nmemory in MB total = %e  global = %e\n\n',tbytes/1e6,gbytes/1e6));

[ss,sr]=system('ps v | grep MATLAB');
printlog(sprintf('\n system memory (xfemfunc): %s\n\n',sr));

