function WriteExoMatrix(exofile,matrix)
%
% write sparse matrix to exodus
%

% % determine maximum bandwidth of matrix
% bmax=0;
% for ir=1:size(matrix,1)
%      [ia,~]=find(matrix(ir,:));
%      bmax=max(length(ia),bmax);
% end
% 
% numrow = size(matrix,1);
% numcol = bmax;
% 
% ncsize=zeros(numrow,1);
% cindex=zeros(numrow,numcol);
% values=zeros(numrow,numcol);
% 
% for ir=1:size(matrix,1)
%      [ia,ib,iv]=find(matrix(ir,:));
%      nc=length(ia);
%      ncsize(ir)=nc;
%      cindex(ir,1:nc)=ib;
%      values(ir,1:nc)=iv;
% end
% 
% % write netcdf file
% 
% ncid = netcdf.create(exofile,'NC_WRITE');
% 
% dimidrow = netcdf.defDim(ncid,'rows',   numrow);
% dimidcol = netcdf.defDim(ncid,'columns',numcol);
% dimidone = netcdf.defDim(ncid,'singlec',1);
% 
% varid_siz = netcdf.defVar(ncid,'ncsize','NC_DOUBLE',[dimidrow dimidone]);
% varid_ind = netcdf.defVar(ncid,'cindex','NC_DOUBLE',[dimidrow dimidcol]);
% varid_val = netcdf.defVar(ncid,'values','NC_DOUBLE',[dimidrow dimidcol]);
% 
% netcdf.endDef(ncid);
% 
% netcdf.putVar(ncid,varid_siz,ncsize);
% netcdf.putVar(ncid,varid_ind,cindex);
% netcdf.putVar(ncid,varid_val,values);

[ia,ib,iv] = find(matrix);

sparseMatVals = [ia,ib,iv];


ncid = netcdf.create(exofile,'NC_WRITE');

dimidtwo = netcdf.defDim(ncid,'matdimslength',2);
dimidone = netcdf.defDim(ncid,'singlec',1);

dimidrow = netcdf.defDim(ncid,'rows',size(sparseMatVals,1));
dimidcol = netcdf.defDim(ncid,'columns',size(sparseMatVals,2));

varid_val = netcdf.defVar(ncid,'sparsemat','NC_DOUBLE',[dimidrow dimidcol]);
varid_fullmatsize = netcdf.defVar(ncid,'fullmatsize','NC_INT',[dimidtwo,dimidone]);

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid_val,double(sparseMatVals));
netcdf.putVar(ncid,varid_fullmatsize,[size(matrix,1),size(matrix,2)]);

netcdf.close(ncid);
end