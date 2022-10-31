function printlog(msgstr)

fid = fopen('matlab.log','a');

if fid < 0
    t = 0;
    while fid < 0 
        t = t+1;
        pause(0.25)
        fid = fopen('matlab.log','a');
        if t > 10
            error('opening log-file failed')
        end
    end
end

fprintf(fid,msgstr);
fclose(fid);
