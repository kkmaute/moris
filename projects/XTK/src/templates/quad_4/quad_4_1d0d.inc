//quad_4_1d0d.inc

moris::Mat<moris::uint> tIndices1d0d(1,8,0);
moris::Mat<moris::uint> tOffsets1d0d(1,5,0);

//Edge 1
tIndices1d0d(0,0) = mNodesInd(0,0);  tIndices1d0d(0,1) = mNodesInd(1,0);
tOffsets1d0d(0,1) = 2;
//Edge 2
tIndices1d0d(0,2) = mNodesInd(1,0);  tIndices1d0d(0,3) = mNodesInd(2,0);
tOffsets1d0d(0,2) = 4;
//Edge 3
tIndices1d0d(0,4) = mNodesInd(2,0);  tIndices1d0d(0,5) = mNodesInd(3,0);
tOffsets1d0d(0,3) = 6;
//Edge 4
tIndices1d0d(0,6) = mNodesInd(0,0);  tIndices1d0d(0,7) = mNodesInd(3,0);
tOffsets1d0d(0,4) = 8;
