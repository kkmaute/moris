    // 1D to 0D ---------------------------------------------------------------------------------------------------------------

    moris::Mat<moris::uint> tIndices1d0d(1,100);
    moris::Mat<moris::uint> tOffsets1d0d(1,51);
    tIndices1d0d(0,0)   =  0 ;    tIndices1d0d(0,1)   =  1 ;// Edge 1
    tIndices1d0d(0,2)   =  1 ;    tIndices1d0d(0,3)   =  2 ;// Edge 2
    tIndices1d0d(0,4)   =  3 ;    tIndices1d0d(0,5)   =  2 ;// Edge 3
    tIndices1d0d(0,6)   =  0 ;    tIndices1d0d(0,7)   =  3 ;// Edge 4
    tIndices1d0d(0,8)   =  4 ;    tIndices1d0d(0,9)   =  5 ;// Edge 5
    tIndices1d0d(0,10)  =  5 ;    tIndices1d0d(0,11)  =  6 ;// Edge 6
    tIndices1d0d(0,12)  =  7 ;    tIndices1d0d(0,13)  =  6 ;// Edge 7
    tIndices1d0d(0,14)  =  4 ;    tIndices1d0d(0,15)  =  7 ;// Edge 8
    tIndices1d0d(0,16)  =  0 ;    tIndices1d0d(0,17)  =  4 ;// Edge 9
    tIndices1d0d(0,18)  =  1 ;    tIndices1d0d(0,19)  =  5 ;// Edge 10
    tIndices1d0d(0,20)  =  2 ;    tIndices1d0d(0,21)  =  6 ;// Edge 11
    tIndices1d0d(0,22)  =  3 ;    tIndices1d0d(0,23)  =  7 ;// Edge 12
    tIndices1d0d(0,24)  =  3 ;    tIndices1d0d(0,25)  =  14;// Edge 13
    tIndices1d0d(0,26)  =  7 ;    tIndices1d0d(0,27)  =  14;// Edge 14
    tIndices1d0d(0,28)  =  6 ;    tIndices1d0d(0,29)  =  14;// Edge 15
    tIndices1d0d(0,30)  =  2 ;    tIndices1d0d(0,31)  =  14;// Edge 16
    tIndices1d0d(0,32)  =  5 ;    tIndices1d0d(0,33)  =  14;// Edge 17
    tIndices1d0d(0,34)  =  1 ;    tIndices1d0d(0,35)  =  14;// Edge 18
    tIndices1d0d(0,36)  =  0 ;    tIndices1d0d(0,37)  =  14;// Edge 19
    tIndices1d0d(0,38)  =  4 ;    tIndices1d0d(0,39)  =  14;// Edge 20
    tIndices1d0d(0,40)  =  3 ;    tIndices1d0d(0,41)  =  10;// Edge 21
    tIndices1d0d(0,42)  =  7 ;    tIndices1d0d(0,43)  =  10;// Edge 22
    tIndices1d0d(0,44)  =  6 ;    tIndices1d0d(0,45)  =  10;// Edge 23
    tIndices1d0d(0,46)  =  2 ;    tIndices1d0d(0,47)  =  10;// Edge 24
    tIndices1d0d(0,48)  =  6 ;    tIndices1d0d(0,49)  =  9 ;// Edge 25
    tIndices1d0d(0,50)  =  2 ;    tIndices1d0d(0,51)  =  9 ;// Edge 26
    tIndices1d0d(0,52)  =  1 ;    tIndices1d0d(0,53)  =  9 ;// Edge 27
    tIndices1d0d(0,54)  =  5 ;    tIndices1d0d(0,55)  =  9 ;// Edge 28
    tIndices1d0d(0,56)  =  4 ;    tIndices1d0d(0,57)  =  8 ;// Edge 29
    tIndices1d0d(0,58)  =  5 ;    tIndices1d0d(0,59)  =  8 ;// Edge 30
    tIndices1d0d(0,60)  =  1 ;    tIndices1d0d(0,61)  =  8 ;// Edge 31
    tIndices1d0d(0,62)  =  0 ;    tIndices1d0d(0,63)  =  8 ;// Edge 32
    tIndices1d0d(0,64)  =  7 ;    tIndices1d0d(0,65)  =  11;// Edge 33
    tIndices1d0d(0,66)  =  3 ;    tIndices1d0d(0,67)  =  11;// Edge 34
    tIndices1d0d(0,68)  =  0 ;    tIndices1d0d(0,69)  =  11;// Edge 35
    tIndices1d0d(0,70)  =  4 ;    tIndices1d0d(0,71)  =  11;// Edge 36
    tIndices1d0d(0,72)  =  7 ;    tIndices1d0d(0,73)  =  13;// Edge 37
    tIndices1d0d(0,74)  =  6 ;    tIndices1d0d(0,75)  =  13;// Edge 38
    tIndices1d0d(0,76)  =  5 ;    tIndices1d0d(0,77)  =  13;// Edge 39
    tIndices1d0d(0,78)  =  4 ;    tIndices1d0d(0,79)  =  13;// Edge 40
    tIndices1d0d(0,80)  =  3 ;    tIndices1d0d(0,81)  =  12;// Edge 41
    tIndices1d0d(0,82)  =  2 ;    tIndices1d0d(0,83)  =  12;// Edge 42
    tIndices1d0d(0,84)  =  1 ;    tIndices1d0d(0,85)  =  12;// Edge 43
    tIndices1d0d(0,86)  =  0 ;    tIndices1d0d(0,87)  =  12;// Edge 44
    tIndices1d0d(0,88)  =  10;    tIndices1d0d(0,89)  =  14;// Edge 45
    tIndices1d0d(0,90)  =  9 ;    tIndices1d0d(0,91)  =  14;// Edge 46
    tIndices1d0d(0,92)  =  8 ;    tIndices1d0d(0,93)  =  14;// Edge 47
    tIndices1d0d(0,94)  =  11;    tIndices1d0d(0,95)  =  14;// Edge 48
    tIndices1d0d(0,96)  =  13;    tIndices1d0d(0,97)  =  14;// Edge 49
    tIndices1d0d(0,98)  =  12;    tIndices1d0d(0,99)  =  14;// Edge 50

    //Offsets
    tOffsets1d0d(0,0)   = 0  ;    tOffsets1d0d(0,1)   = 2  ;    tOffsets1d0d(0,2)   = 4  ;    tOffsets1d0d(0,3)   = 6;
    tOffsets1d0d(0,4)   = 8  ;    tOffsets1d0d(0,5)   = 10 ;    tOffsets1d0d(0,6)   = 12 ;    tOffsets1d0d(0,7)   = 14;
    tOffsets1d0d(0,8)   = 16 ;    tOffsets1d0d(0,9)   = 18 ;    tOffsets1d0d(0,10)  = 20 ;    tOffsets1d0d(0,11)  = 22;
    tOffsets1d0d(0,12)  = 24 ;    tOffsets1d0d(0,13)  = 26 ;    tOffsets1d0d(0,14)  = 28 ;    tOffsets1d0d(0,15)  = 30;
    tOffsets1d0d(0,16)  = 32 ;    tOffsets1d0d(0,17)  = 34 ;    tOffsets1d0d(0,18)  = 36 ;    tOffsets1d0d(0,19)  = 38;
    tOffsets1d0d(0,20)  = 40 ;    tOffsets1d0d(0,21)  = 42 ;    tOffsets1d0d(0,22)  = 44 ;    tOffsets1d0d(0,23)  = 46;
    tOffsets1d0d(0,24)  = 48 ;    tOffsets1d0d(0,25)  = 50 ;    tOffsets1d0d(0,26)  = 52 ;    tOffsets1d0d(0,27)  = 54;
    tOffsets1d0d(0,28)  = 56 ;    tOffsets1d0d(0,29)  = 58 ;    tOffsets1d0d(0,30)  = 60 ;    tOffsets1d0d(0,31)  = 62;
    tOffsets1d0d(0,32)  = 64 ;    tOffsets1d0d(0,33)  = 66 ;    tOffsets1d0d(0,34)  = 68 ;    tOffsets1d0d(0,35)  = 70;
    tOffsets1d0d(0,36)  = 72 ;    tOffsets1d0d(0,37)  = 74 ;    tOffsets1d0d(0,38)  = 76 ;    tOffsets1d0d(0,39)  = 78;
    tOffsets1d0d(0,40)  = 80 ;    tOffsets1d0d(0,41)  = 82 ;    tOffsets1d0d(0,42)  = 84 ;    tOffsets1d0d(0,43)  = 86;
    tOffsets1d0d(0,44)  = 88 ;    tOffsets1d0d(0,45)  = 90 ;    tOffsets1d0d(0,46)  = 92 ;    tOffsets1d0d(0,47)  = 94;
    tOffsets1d0d(0,48)  = 96 ;    tOffsets1d0d(0,49)  = 98 ;    tOffsets1d0d(0,50)  = 100;
