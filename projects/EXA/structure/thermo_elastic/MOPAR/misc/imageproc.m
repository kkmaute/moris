function [levelset,ex,ey,nx,ny]=imageproc(imgfil,nx,ny,plx,cxy,invflag,pltflag,offset,nsmooth,ninit)
%function [levelset,ex,ey,nx,ny]=imageproc(imgfil,nx,ny,plx,cxy,invflag,pltflag,offset,nsmooth,ninit)
%
% create level set represenation from image file
% 
% input:  imgfil  -  name of impage file (prefered: png format)
%         nx,ny   -  size of levelset mesh
%         plx     -  physical size of image (only x-direction can be specified)
%         cxy     -  cxy=[cx,cy]; number of pixels cut off from boundaries 
%                    to avoid fringing effects (def: [0,0])
%         invflag -  invert image 0: no(def), 1: yes
%         pltflag -  print figures during image processing (def: no)
%         offset  -  offset of level set function [0;1] (def: 0.5) 
%         nsmooth -  smoothing level (def: 50)
%         ninit   -  iterations in creating signed distance function (def: 50) 
%
% one can check results with:
%
%  contour(ex,ey,reshape(levelset,ny,nx));hold on;
%  contour(ex,ey,reshape(levelset,ny,nx),[0 0],'LineWidth',2)
%
addpath(genpath(pwd));

if nargin < 10; ninit   = 500;  end
if nargin <  9; nsmooth = 50;  end
if nargin <  8; offset  = 0.5; end
if nargin <  7; pltflag = 0;   end
if nargin <  6; invflag = 0;   end
if nargin <  5; cxy=[1 1];     end

% read image
im = imread(imgfil);
phi = double(im(:,:,1));

% crop file by cx and cy
cx=cxy(1); cy=cxy(2);
[lxp lyp]=size(phi);
phi=phi(cx:lxp-2*cx,cy:lyp-2*cy);

% invert figure
if invflag; phi=-phi;end

if pltflag; figure; imagesc(phi); axis image; colormap gray; end
if pltflag; figure; contour(phi); axis image; end

% b:weighting for curvature-based force. b needs to be positive
b = 0.3*ones(size(phi));

dx=1;
dy=1;
phi = evolve2D(phi,dx,dy,0.5,nsmooth,[],[],0,[],0,[],[],1,b);

if pltflag; figure; imagesc(phi); axis image; colormap gray; end
if pltflag; figure; contour(phi); axis image; end

% scale level set function to be [0,1]
% shift level set function such that zero-level shows grain boundaries

minp=min(min(phi));
maxp=max(max(phi));
phiscl=1/(maxp-minp)*(phi-minp*ones(size(phi)));
phiscl=(phiscl-offset*ones(size(phi)));

if pltflag; figure; contour(phiscl,[0 0],'b','linewidth',1); axis image; end

% create signed distace approximation

phiscls=reinit_SD(phiscl,dx,dy,0.8,'ENO1',ninit);

% interpolate on given mesh (watch node numbering)

ly=size(phi,1);
lx=size(phi,2);

dl=lx/nx;
if isempty(ny)
    ny=ceil(ly/dl);
end

[x,y]   = meshgrid(1:lx,1:ly);
[ex,ey] = meshgrid(1:nx,1:ny);
ex      = dl*ex;
ey      = dl*ey;

phie = interp2(x,y,phiscls,ex,ey);

for i=1:nx
    for j=1:ny
        if isnan(phie(j,i)); phie(j,i)=0;end
    end
end

if pltflag; figure; contour(ex,ey,phie,[0 0],'b','linewidth',1); axis image; end
if pltflag; figure; meshc(x,y,min(lx,ly)/max(max(phiscls))*phiscls); axis image; end

% write file

sclx=plx/dl/max(nx,ny);

fid=fopen('levelset.data','w');

fprintf(fid,'# node-id  x-coord  y-coord  levelset\n');
in=0;
for ix=1:ny
    for iy=1:nx
        in=in+1;
        fprintf(fid,'%d  %f  %f  %f\n',in,sclx*(ex(ix,iy)-1),sclx*(ey(ix,iy)-1),phie(ix,iy));
    end
end
fclose(fid);

if pltflag; 
    figure; imagesc(phi); axis image; colormap gray;
    hold on;
    contour(x,y,phiscls,[0 0],'b','linewidth',1); axis image;
    contour(ex,ey,phie,[0 0],'r','linewidth',1); axis image;
end

% shift and scale mesh before output

ex=ex-1;
ey=ey-1;

ex=sclx*ex;
ey=sclx*ey;

levelset=reshape(phie,nx*ny,1);
save levelset.mat levelset nx ny ex ey
