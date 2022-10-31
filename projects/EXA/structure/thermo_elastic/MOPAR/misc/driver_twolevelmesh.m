function driver_twolevelmesh

xlen=20;   % dimension of design domain in x-direction
ylen=10;   % dimension of design domain in y-direction

nx_fine=20;  % fine mesh: intervalls (x-direction)
ny_fine=10;  % fine mesh: intervalls (y-direction)

nx_crse=5;   % coarse mesh: intervalls (x-direction)
ny_crse=2;   % coarse mesh: intervalls (y-direction)

[fecoord_fine,fetopo_fine,fecoord_crse,fetopo_crse,Mmat,bmat]= ...
    twolevelmesh(xlen,ylen,nx_fine,ny_fine,nx_crse,ny_crse);

phi_crse=fecoord_crse(:,1);
phi_fine=Mmat\(bmat*phi_crse);

VarName{1} = 'phi';

CreateNetcdfFile('coarse.exo',{'quad4'},{'x','y'},fecoord_crse,fetopo_crse,size(fetopo_crse,1),zeros(size(fecoord_crse,1),1),1,size(fetopo_crse,1));
IniExoData('coarse.exo',fecoord_crse);
WriteExoData('coarse.exo',VarName,phi_crse,1,0.0);


CreateNetcdfFile('fine.exo',{'quad4'},{'x','y'},fecoord_fine,fetopo_fine,size(fetopo_fine,1),zeros(size(fecoord_fine,1),1),1,size(fetopo_fine,1));
IniExoData('fine.exo',fecoord_fine);
WriteExoData('fine.exo',VarName,phi_fine,1,0.0);

end
