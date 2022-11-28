function savecell2file(data,dname,fname,sflag)

% save data to file in cell format
%
% sflag - 0 save to new file
%       - 1 append to file
%

if sflag==1
    if exist(fname)
        fnum=length(fieldnames(load(fname)));
        for n=1:size(data,2)
            d=dname{n};
            eval([d '_' num2str(fnum+1) '=data{n};']);
            save(fname,'-append','-regexp',[d '_*'])
        end
    else
        for n=1:size(data,2)
            d=dname{n};
            eval([d '_' num2str(1) '=data{n};']);
            if n==1
                save(fname,'-regexp',[d '_*'])
            else
                save(fname,'-append','-regexp',[d '_*'])
            end
        end
    end
elseif sflag==0
    for n=1:size(data,2)
        d=dname{n};
        eval([d '_' num2str(1) '=data{n};']);
        if n==1
            save(fname,'-regexp',[d '_*'])
        else
            save(fname,'-append','-regexp',[d '_*'])
        end
    end
end

