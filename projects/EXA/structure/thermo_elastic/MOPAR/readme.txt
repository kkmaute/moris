* edit techplot file and change name of variables to x and y
* load techplot file in paraview and write a exodus files (both materials): allmat.e
* load both materials exodus file into paraview 
  - threshold by matNum [2,2] and write mat2.e
  - invert thresold and write mat1.e
* load mat1.e and mat2.e into cubit

*------------------------------------------------------

import mesh geometry "/data/home/maute/work/Sam/MOPAR/2021-11-18/45km_2mat_abl_new/mat1.e" feature_angle 135.00 merge 
import mesh geometry "/data/home/maute/work/Sam/MOPAR/2021-11-18/45km_2mat_abl_new/mat2.e" feature_angle 135.00 merge 
set duplicate block elements off
reset block
set duplicate block elements off
block 1 add surface 5  
set duplicate block elements off
block 2 add surface 1  
block 1 2 element type tri3
Sideset 1 add curve 2  
Sideset 2 add curve 3 5 
block 1 name "block 1"
block 2 name "block 2"
Sideset 1 name "side 1"
Sideset 2 name "side 2"
set exodus netcdf4 off
set large exodus on
export mesh "/data/home/maute/work/Sam/MOPAR/2021-11-18/45km_2mat_abl_new/exomesh.e" dimension 2 overwrite
reset
import mesh geometry "/data/home/maute/work/Sam/MOPAR/2021-11-18/45km_2mat_abl_new/exomesh.e" feature_angle 135.00 merge merge_nodes  
topology check coincident node surface 2 1  tolerance 1.0e-6 draw brief result group 
set exodus netcdf4 off
set large exodus on
export mesh "/data/home/maute/work/Sam/MOPAR/2021-11-18/45km_2mat_abl_new/exomesh.e" dimension 2 overwrite

*------------------------------------------------------

* rm history* cubit*   
* pwd
* cp exomesh.e ../Mapper/.
* cp allmat.e ../Mapper/.
* cd ../Mapper/
* /home/MATLAB/R2020b/bin/matlab
* cd back to working directory (see above)
* cp ../Mapper/exomesh.e .
* moris.sh opt 1 MOPAR_Problem

