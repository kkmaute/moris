// 0D to 1D ---------------------------------------------------------------------------------------------------------------
    moris::Mat<moris::uint> tIndices0d1d(1,94);
    moris::Mat<moris::uint> tOffsets0d1d(1,16);

    //Node 1 (n1)
    tIndices0d1d(0,0)  = 1 ;    tIndices0d1d(0,1)  = 4 ;   tIndices0d1d(0,2)  = 9 ;  tIndices0d1d(0,3)  = 19;
    tIndices0d1d(0,4)  = 32;    tIndices0d1d(0,5)  = 35;   tIndices0d1d(0,6)  = 44;

    //Node 2 (n2)
    tIndices0d1d(0,7)  = 1 ;    tIndices0d1d(0,8)  = 2 ;   tIndices0d1d(0,9)  = 10;  tIndices0d1d(0,10) = 18;
    tIndices0d1d(0,11) = 27;    tIndices0d1d(0,12) = 31;   tIndices0d1d(0,13) = 43;

    //Node 3 (n3)
    tIndices0d1d(0,14) = 1 ;    tIndices0d1d(0,15) = 2 ;   tIndices0d1d(0,16) = 10;  tIndices0d1d(0,17) = 18;
    tIndices0d1d(0,18) = 27;    tIndices0d1d(0,19) = 31;   tIndices0d1d(0,20) = 43;

    //Node 4 (n4)
    tIndices0d1d(0,21) = 2 ;    tIndices0d1d(0,22) = 3 ;    tIndices0d1d(0,23) = 11;    tIndices0d1d(0,24) = 16;
    tIndices0d1d(0,25) = 24;    tIndices0d1d(0,26) = 26;    tIndices0d1d(0,27) = 42;

    //Node 5 (n5)
    tIndices0d1d(0,28) = 5 ;    tIndices0d1d(0,29) = 8 ;    tIndices0d1d(0,30) = 9 ;    tIndices0d1d(0,31) = 20;
    tIndices0d1d(0,32) = 29;    tIndices0d1d(0,33) = 36;    tIndices0d1d(0,34) = 40;

    //Node 6 (n6)
    tIndices0d1d(0,35) = 5 ;    tIndices0d1d(0,36) = 6 ;    tIndices0d1d(0,37) = 10 ;   tIndices0d1d(0,38) = 17;
    tIndices0d1d(0,39) = 28;    tIndices0d1d(0,40) = 30;    tIndices0d1d(0,41) = 39;

    //Node 7 (n7)
    tIndices0d1d(0,42) = 6 ;    tIndices0d1d(0,43) = 7 ;    tIndices0d1d(0,44) = 11 ;   tIndices0d1d(0,45) = 15;
    tIndices0d1d(0,46) = 23;    tIndices0d1d(0,47) = 25;    tIndices0d1d(0,48) = 38;

    //Node 8 (n8)
    tIndices0d1d(0,49) = 7 ;    tIndices0d1d(0,50) = 8 ;    tIndices0d1d(0,51) = 12 ;   tIndices0d1d(0,52) = 14;
    tIndices0d1d(0,53) = 22;    tIndices0d1d(0,54) = 33;    tIndices0d1d(0,55) = 37;

    //Node 9 (n9)
    tIndices0d1d(0,56) = 29;    tIndices0d1d(0,57) = 30;    tIndices0d1d(0,58) = 31;    tIndices0d1d(0,59) = 32;

    //Node 10 (n10)
    tIndices0d1d(0,60) = 25;    tIndices0d1d(0,61) = 26;    tIndices0d1d(0,62) = 27;    tIndices0d1d(0,63) = 28;

    //Node 11 (n11)
    tIndices0d1d(0,64) = 21;    tIndices0d1d(0,65) = 22;    tIndices0d1d(0,66) = 23;    tIndices0d1d(0,67) = 24;

    //Node 12 (n12)
    tIndices0d1d(0,68) = 33;    tIndices0d1d(0,69) = 34;    tIndices0d1d(0,70) = 35;    tIndices0d1d(0,71) = 36;

    //Node 13 (n13)
    tIndices0d1d(0,72) = 41;    tIndices0d1d(0,73) = 42;    tIndices0d1d(0,74) = 43;    tIndices0d1d(0,75) = 44;

    //Node 14 (n14)
    tIndices0d1d(0,76) = 37;    tIndices0d1d(0,77) = 38;    tIndices0d1d(0,78) = 39;    tIndices0d1d(0,79) = 40;

    //Node 15 (n15)
    tIndices0d1d(0,80) = 13;    tIndices0d1d(0,81) = 14;    tIndices0d1d(0,82) = 15;    tIndices0d1d(0,83) = 16;
    tIndices0d1d(0,84) = 17;    tIndices0d1d(0,85) = 18;    tIndices0d1d(0,86) = 19;    tIndices0d1d(0,87) = 20;
    tIndices0d1d(0,88) = 45;    tIndices0d1d(0,89) = 46;    tIndices0d1d(0,90) = 47;    tIndices0d1d(0,91) = 48;
    tIndices0d1d(0,92) = 49;    tIndices0d1d(0,93) = 50;


    //Offsets
    tOffsets0d1d(0,0)  = 0 ;    tOffsets0d1d(0,1)  = 7 ;    tOffsets0d1d(0,2)  = 14;    tOffsets0d1d(0,3) = 21;
    tOffsets0d1d(0,4)  = 28;    tOffsets0d1d(0,5)  = 35;    tOffsets0d1d(0,6)  = 42;    tOffsets0d1d(0,7) = 49;
    tOffsets0d1d(0,8)  = 56;    tOffsets0d1d(0,9)  = 60;    tOffsets0d1d(0,10) = 64;    tOffsets0d1d(0,11) = 68;
    tOffsets0d1d(0,12) = 72;    tOffsets0d1d(0,13) = 76;    tOffsets0d1d(0,14) = 80;    tOffsets0d1d(0,15) = 94;
