* edit techplot file and change name of variables to x and y
* load techplot file in paraview and write a exodus files (both materials): allmat.e
* load both materials exodus file into paraview 
  - threshold by matNum [2,2] and write mat2.e
  - invert thresold and write mat1.e
* load mat1.e and mat2.e into cubit

*------------------------------------------------------

import mesh geometry "mat1.e" feature_angle 135.00 merge 
import mesh geometry "mat2.e" feature_angle 135.00 merge 
set duplicate block elements off
reset block
set duplicate block elements off
block 1 add surface 5  
set duplicate block elements off
block 2 add surface 1  
block 1 2 element type tri3
Sideset 1 add curve 2  
Sideset 2 add curve 3 
Sideset 3 add curve 5 
block 1 name "block 1"
block 2 name "block 2"
Sideset 1 name "side 1"
Sideset 2 name "side 2"
Sideset 3 name "side 3"
set exodus netcdf4 off
set large exodus on
export mesh "exomesh.e" dimension 2 overwrite
reset
import mesh geometry "exomesh.e" feature_angle 135.00 merge merge_nodes  
topology check coincident node surface 3 2 1  tolerance 1.0e-6 draw brief result group 
set exodus netcdf4 off
set large exodus on
export mesh "exomesh.e" dimension 2 overwrite

*------------------------------------------------------

* rm history* cubit*   

* run matlab script: mapper.m to generate 

* cp exomesh.e traction.e

