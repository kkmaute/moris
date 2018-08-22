//quad_4_1d2d.inc

moris::Mat<moris::uint> tIndices1d2d(1,4,0);
moris::Mat<moris::uint> tOffsets1d2d(1,5,0);

//Edge 1
tIndices1d2d(0,0) = 1;
tOffsets1d2d(0,1) = 2;
//Edge 2
tIndices1d2d(0,1) = 1;
tOffsets1d2d(0,2) = 3;
//Edge 3
tIndices1d2d(0,2) = 1;
tOffsets1d2d(0,3) = 4;
//Edge 4
tIndices1d2d(0,3) = 1;
tOffsets1d2d(0,4) = 4;


