%
%   Copyright (c) 2023 University of Colorado
%   Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
%
%  ------------------------------------------------------------------------------------
%
function Image_to_Levelset_3D
%==========================================================================
% script to convert stack of bitmap images into 3D signed-distance level
% set file for moris image_sdf_field input
%
% input (set manually below)
%
%        1. input image file name (e.g. filename.xxx.png)
%        2. output file name (hdf5, without suffix)
%        3. rotate angle in x-y plane
%        4. padding size in voxels
%        5. blurring factor (smoothes picture)
%
% output
%
%        1. signed distance level set field (hdf5 format)
%        2. png files with signed distance field (can be read in paraview)
%           Img3D_tmp.<id number>.png
%
%==========================================================================

close all
clear all

% image file name
fname="input_image_file.png";

% background mesh file name (without file ending)
outputbase="output_file_name";

% rotation angle in degrees
rangle=0;

% buffer size with user-defined or additional phase
bsize=0;

% smoothing factor
blur=0;

% number of processors moris will be executed on
nproc=1;

% scaling factor for the resulting image SDF
SDF_scale = 1.0;

%==========================================================================
% no need to edit code below
clc;

delete(sprintf("%s_*.hdf5",outputbase));

% process image
image=-SDF_scale*procimage(fname,rangle,blur,bsize);

% save to hdf5
savetohdf5(image,outputbase,nproc);
end

%==========================================================================

function savetohdf5(image,outputbase,nproc)

% get dimensions
[nx,ny,nz]=size(image);

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
    h5create(outfile, dataset, [3 1]);
    h5write(outfile, dataset, [nx; ny; nz]);

    dataset=sprintf('%s%s','/','SDF');
    h5create(outfile, dataset, [nx*ny*nz 1]);
    h5write(outfile, dataset, imgvec);
end
end

%==========================================================================

function [image3D]=procimage(fname,rangle,blur,bsize)

% determine number of files
filelist = dir(sprintf("%s.*.*",fname));

% get number of files
nz=size(filelist,1);

for iz=1:nz

    % name of 2D image
    fname=filelist(iz).name;

    % read image
    im=imread(fname);

    % rotate image
    im=imrotate(im,rangle);

    % get size of 2D image
    [nx,ny]=size(im);

    % convert to bw image
    if size(im,3) > 1
        image=rgb2gray(im);
        imageinv=255-image;
    else
        imageinv=255-im;
    end

    % allocate 3D image array
    if iz==1
        fprintf("\n\nVoxel dimension (original) nx=%d  ny=%d  nz=%d  x/y-ratio=%f  x/z-ratio=%f\n", ...
            nx,ny,nz,nx/ny,nx/nz);
        
        bvalue=double(imageinv(1,1));
        image3D=uint8(bvalue*ones(nx+2*bsize,ny+2*bsize,nz+2*bsize));

        fprintf("\nVoxel dimension (original) nx=%d  ny=%d  nz=%d  x/y-ratio=%f  x/z-ratio=%f\n\n", ...
            nx+2*bsize,ny+2*bsize,nz+2*bsize,(nx+2*bsize)/(ny+2*bsize),(nx+2*bsize)/(nz+2*bsize));
    end

    % store 2D image in 3D array
    image3D(bsize+1:bsize+nx,bsize+1:bsize+ny,bsize+iz)=imageinv;

end

% smooth image
if (blur > 0)
    image3D=imgaussfilt3(image3D,blur);
end

% threshold
imgtmp=double(image3D)/255;
imgthr=double(imgtmp>0.5);

% compute distance field
image3D=bwdist(imgthr)-bwdist(1-imgthr);

isosurface(image3D,0), axis equal, view(3)
camlight, lighting gouraud

% output stack of images with signed distance field
delete Img3D_tmp.*.png

minval=min(min(min(image3D)));
maxval=max(max(max(image3D)));

for iz=1:nz+2*bsize
    fname=sprintf("Img3D_tmp.%03d.png",iz);
    imwrite(250*(image3D(:,:,iz)-minval)/(maxval-minval),colormap("jet"),fname)
end

end
