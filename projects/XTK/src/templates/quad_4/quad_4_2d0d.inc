//quad_4_2d0d.inc


moris::Mat<moris::uint> tIndices2d0d(1,4);
moris::Mat<moris::uint> tOffsets2d0d(1,2);

// element 1
tIndices2d0d(0,0) = mNodesInd(0,0);    tIndices2d0d(0,1) = mNodesInd(1,0);    tIndices2d0d(0,2) = mNodesInd(2,0);    tIndices2d0d(0,3) = mNodesInd(3,0);
tOffsets2d0d(0,0) = 0;    tOffsets2d0d(0,1) = 4;
