    moris::Mat<moris::uint> tIndices0d1d(1,8,0);
    moris::Mat<moris::uint> tOffsets0d1d(1,5,0);

    // Node 1
    tIndices0d1d(0,0) = 1;    tIndices0d1d(0,1) = 4;
    tOffsets0d1d(0,1) = 2;
    // Node 2
    tIndices0d1d(0,2) = 1;    tIndices0d1d(0,3) = 2;
    tOffsets0d1d(0,2) = 4;
    // Node 3
    tIndices0d1d(0,4) = 2;    tIndices0d1d(0,5) = 3;
    tOffsets0d1d(0,3) = 6;
    // Node 4
    tIndices0d1d(0,6) = 3;    tIndices0d1d(0,7) = 4;
    tOffsets0d1d(0,4) = 8;
