moris::Mat<moris::uint> tAnc2D(1,60,UINT_MAX);

// XTK faces created on parent element local face 1
tAnc2D(0,8)  = aParentEntities(2)(0,0);
tAnc2D(0,9)  = aParentEntities(2)(0,0);
tAnc2D(0,10) = aParentEntities(2)(0,0);
tAnc2D(0,11) = aParentEntities(2)(0,0);

// XTK faces created on parent element local face 2
tAnc2D(0,4) = aParentEntities(2)(0,1);
tAnc2D(0,5) = aParentEntities(2)(0,1);
tAnc2D(0,6) = aParentEntities(2)(0,1);
tAnc2D(0,7) = aParentEntities(2)(0,1);

// XTK faces created on parent element local face 3
tAnc2D(0,0) = aParentEntities(2)(0,2);
tAnc2D(0,1) = aParentEntities(2)(0,2);
tAnc2D(0,2) = aParentEntities(2)(0,2);
tAnc2D(0,3) = aParentEntities(2)(0,2);

// XTK faces created on parent element local face 4
tAnc2D(0,12) = aParentEntities(2)(0,3);
tAnc2D(0,13) = aParentEntities(2)(0,3);
tAnc2D(0,14) = aParentEntities(2)(0,3);
tAnc2D(0,15) = aParentEntities(2)(0,3);

// XTK faces created on parent element local face 5
tAnc2D(0,20) = aParentEntities(2)(0,4);
tAnc2D(0,21) = aParentEntities(2)(0,4);
tAnc2D(0,22) = aParentEntities(2)(0,4);
tAnc2D(0,23) = aParentEntities(2)(0,4);

// XTK faces created on parent element local face 6
tAnc2D(0,16) = aParentEntities(2)(0,5);
tAnc2D(0,17) = aParentEntities(2)(0,5);
tAnc2D(0,18) = aParentEntities(2)(0,5);
tAnc2D(0,19) = aParentEntities(2)(0,5);

// Everything else has UINT_MAX indicating it does not have a parent face entity
