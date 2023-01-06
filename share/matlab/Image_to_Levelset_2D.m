% 
%   Copyright (c) 2023 University of Colorado
%   Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%  
%  ------------------------------------------------------------------------------------
%  
function Image_to_Levelset_2D
%==========================================================================
% script to convert bitmap image into signed-distance level set file for
% moris image_sdf_field input
%
% input (set manually below)
%
%        1. input image file name
%        2. output file name (hdf5, without suffix)
%        3. rotate angle 
%        4. padding size in voxels
%        5. blurring factor (smoothes picture)
%
% output
%
%        1. signed distance level set field (hdf5 format)
%==========================================================================

close all
clear all

% image file name
fname="oscillatorForKurt2_Trimmed.bmp";

% background mesh file name (base)
outputbase="Oscillator_real.HMR";

% rotation angle in degrees
rangle=-90;

% buffer size with user-defined or additional phase
bsize=20;

% smoothing factor
blur=0;

% number of processors moris will be executed on
nproc=16;

%==========================================================================
% no need to edit code below
clc;

% process image
image=procimage(fname,rangle,blur,bsize);

% save to hdf5
savetohdf5(image,outputbase,nproc);
end

%==========================================================================

function savetohdf5(image,outputbase,nproc)

% get dimensions
[nx,ny]=size(image);

% reshape image data
imgvec=reshape(image,[],1);

% save data to hdf5
for ip=1:nproc

    if nproc>1
        outfile=sprintf('./%s_%d.%d.hdf5',outputbase,nproc,ip-1);
    else
        outfile=sprintf('./%s.hdf5',outputbase);
    end

    fprintf(' ... writing data to file: %s\n',outfile);

    if exist(outfile,"file")
        delete(outfile)
    end

    dataset=sprintf('%s%s','/','Dimensions');
    h5create(outfile, dataset, [2 1]);
    h5write(outfile, dataset, [nx; ny]);

    dataset=sprintf('%s%s','/','SDF');
    h5create(outfile, dataset, [nx*ny 1]);
    h5write(outfile, dataset, imgvec);
end
end

%==========================================================================

function [image]=procimage(fname,rangle,blur,bsize)

% read image
im=imread(fname);

% rotate image
im=imrotate(im,rangle);

% convert to bw image
if size(im,3) > 1
    image=rgb2gray(im);
    imageinv=255-image;
else
    error('Incorrect image type - needs to be RGB');
end

imshow(imageinv)

[nx,ny]=size(imageinv);
fprintf("\n\n Voxel dimension (original) nx = %d  ny = %d  x/y-ratio = %f\n",nx,ny,nx/ny);

if bsize > 0
    nxb=nx+2*bsize;
    nyb=ny+2*bsize;
    imageinvb=double(imageinv(1,1))*ones(nxb,nyb);
    imageinvb(bsize+1:bsize+nx,bsize+1:bsize+ny)=imageinv;
    imageinv=imageinvb;

    [nx,ny]=size(imageinv);
    fprintf("\n Voxel dimension (with buffer) nx = %d  ny = %d  x/y-ratio = %f\n\n",nx,ny,nx/ny);
end

% smooth image
if (blur > 0)
    image=imgaussfilt(imageinv,blur);
else
    image=imageinv;
end

% threshold 
imgtmp=double(image)/255;
imgthr=double(imgtmp>0.5);

image=bwdist(imgthr)-bwdist(1-imgthr);
contourf(image)
colorbar
axis equal

end
