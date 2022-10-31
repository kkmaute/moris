%
% script to compute error norm of fields
%
% output - original, mapped, and error fields in exodus format
%        - L2 and L_infinity errors
%

clear all;
close all;

exofile='fields.exo';

% data to be mapped 

load MapData.mat

xp=Coords_40x60_XFEM(:,1);
yp=Coords_40x60_XFEM(:,2);
val=pres_40x60_XFEM;

% reference mesh and solution

ref_mesh='Ubend_40x60_bodyfitted.g';
p_org=pres_40x60_bf;

% build interpolation of scattered data

mlv = TriScatteredInterp(xp,yp,val,'linear');

[MassMat,fecoord,fetopo,feprop,febc,numfenod,numfele]=ProcessExoMesh(ref_mesh,0);

p_mapped=mlv(fecoord(:,1),fecoord(:,2));

if sum(find(isnan(p_mapped))) > 0; warning('mapping failed - contains nans'); end 

% compute error and output original and mapped fields

p_error=p_org-p_mapped;

CreateNetcdfFile(exofile,{'quad4'},{'x','y'},fecoord,fetopo,numfele,febc,1,numfele);
IniExoData(exofile,fecoord);

VarName{1} = 'original';
VarName{2} = 'mapped';
VarName{3} = 'local error';

WriteExoData(exofile,VarName,[p_org,p_mapped,abs(p_error)],1,1);

errnorm=(p_error'*MassMat*p_error)/(p_org'*MassMat*p_org)*100;
errmax=max(abs(p_error))/max(p_org)*100;

fprintf('L2 error = %e percent ;   max local error = %e percent \n',errnorm,errmax);
