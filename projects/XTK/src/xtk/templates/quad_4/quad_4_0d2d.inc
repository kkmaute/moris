//quad_4_0d2d.inc

moris::Mat<moris::uint> tIndices0d2d(1,4,0);
moris::Mat<moris::uint> tOffsets0d2d(1,5,0);


// Node 1
tIndices0d2d(0,0) = 1;
tOffsets0d2d(0,1) = 1;
// Node 2
tIndices0d2d(0,1) = 1;
tOffsets0d2d(0,2) = 2;
// Node 3
tIndices0d2d(0,2) = 1;
tOffsets0d2d(0,3) = 3;
// Node 4
tIndices0d2d(0,3) = 1;
tOffsets0d2d(0,4) = 4;
