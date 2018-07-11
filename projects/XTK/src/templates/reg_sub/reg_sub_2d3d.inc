// 2D to 3D ---------------------------------------------------------------------------------------------------------------

    moris::Mat<moris::uint> tIndices2d3d(1,96);
    moris::Mat<moris::uint> tOffsets2d3d(1,61);


    tIndices2d3d(0,0)   = 1 ;                             // Face 1
    tIndices2d3d(0,1)   = 2 ;                             // Face 2
    tIndices2d3d(0,2)   = 3 ;                             // Face 3
    tIndices2d3d(0,3)   = 4 ;                             // Face 4
    tIndices2d3d(0,4)   = 5 ;                             // Face 5
    tIndices2d3d(0,5)   = 6 ;                             // Face 6
    tIndices2d3d(0,6)   = 7 ;                             // Face 7
    tIndices2d3d(0,7)   = 8 ;                             // Face 8
    tIndices2d3d(0,8)   = 9 ;                             // Face 9
    tIndices2d3d(0,9)   = 10;                             // Face 10
    tIndices2d3d(0,10)  = 11;                             // Face 11
    tIndices2d3d(0,11)  = 12;                             // Face 12
    tIndices2d3d(0,12)  = 13;                             // Face 13
    tIndices2d3d(0,13)  = 14;                             // Face 14
    tIndices2d3d(0,14)  = 15;                             // Face 15
    tIndices2d3d(0,15)  = 16;                             // Face 16
    tIndices2d3d(0,16)  = 17;                             // Face 17
    tIndices2d3d(0,17)  = 18;                             // Face 18
    tIndices2d3d(0,18)  = 19;                             // Face 19
    tIndices2d3d(0,19)  = 20;                             // Face 20
    tIndices2d3d(0,20)  = 21;                             // Face 21
    tIndices2d3d(0,21)  = 22;                             // Face 22
    tIndices2d3d(0,22)  = 23;                             // Face 23
    tIndices2d3d(0,23)  = 24;                             // Face 24
    tIndices2d3d(0,24)  = 1 ;    tIndices2d3d(0,25)  = 2 ;// Face 25
    tIndices2d3d(0,26)  = 2 ;    tIndices2d3d(0,27)  = 3 ;// Face 26
    tIndices2d3d(0,28)  = 3 ;    tIndices2d3d(0,29)  = 4 ;// Face 27
    tIndices2d3d(0,30)  = 1 ;    tIndices2d3d(0,31)  = 4 ;// Face 28
    tIndices2d3d(0,32)  = 5 ;    tIndices2d3d(0,33)  = 6 ;// Face 29
    tIndices2d3d(0,34)  = 6 ;    tIndices2d3d(0,35)  = 7 ;// Face 30
    tIndices2d3d(0,36)  = 7 ;    tIndices2d3d(0,37)  = 8 ;// Face 31
    tIndices2d3d(0,38)  = 5 ;    tIndices2d3d(0,39)  = 8 ;// Face 32
    tIndices2d3d(0,40)  = 9 ;    tIndices2d3d(0,41)  = 10;// Face 33
    tIndices2d3d(0,42)  = 10;    tIndices2d3d(0,43)  = 11;// Face 34
    tIndices2d3d(0,44)  = 11;    tIndices2d3d(0,45)  = 12;// Face 35
    tIndices2d3d(0,46)  = 9 ;    tIndices2d3d(0,47)  = 12;// Face 36
    tIndices2d3d(0,48)  = 13;    tIndices2d3d(0,49)  = 14;// Face 37
    tIndices2d3d(0,50)  = 14;    tIndices2d3d(0,51)  = 15;// Face 38
    tIndices2d3d(0,52)  = 15;    tIndices2d3d(0,53)  = 16;// Face 39
    tIndices2d3d(0,54)  = 13;    tIndices2d3d(0,55)  = 16;// Face 40
    tIndices2d3d(0,56)  = 17;    tIndices2d3d(0,57)  = 18;// Face 41
    tIndices2d3d(0,58)  = 18;    tIndices2d3d(0,59)  = 19;// Face 42
    tIndices2d3d(0,60)  = 19;    tIndices2d3d(0,61)  = 20;// Face 43
    tIndices2d3d(0,62)  = 17;    tIndices2d3d(0,63)  = 20;// Face 44
    tIndices2d3d(0,64)  = 21;    tIndices2d3d(0,65)  = 22;// Face 45
    tIndices2d3d(0,66)  = 22;    tIndices2d3d(0,67)  = 23;// Face 46
    tIndices2d3d(0,68)  = 23;    tIndices2d3d(0,69)  = 24;// Face 47
    tIndices2d3d(0,70)  = 21;    tIndices2d3d(0,71)  = 24;// Face 48
    tIndices2d3d(0,72)  = 1 ;    tIndices2d3d(0,73)  = 5 ;// Face 49
    tIndices2d3d(0,74)  = 2 ;    tIndices2d3d(0,75)  = 21;// Face 50
    tIndices2d3d(0,76)  = 3 ;    tIndices2d3d(0,77)  = 13;// Face 51
    tIndices2d3d(0,78)  = 4 ;    tIndices2d3d(0,79)  = 17;// Face 52
    tIndices2d3d(0,80)  = 6 ;    tIndices2d3d(0,81)  = 22;// Face 53
    tIndices2d3d(0,82)  = 7 ;    tIndices2d3d(0,83)  = 9 ;// Face 54
    tIndices2d3d(0,84)  = 8 ;    tIndices2d3d(0,85)  = 18;// Face 55
    tIndices2d3d(0,86)  = 10;    tIndices2d3d(0,87)  = 23;// Face 56
    tIndices2d3d(0,88)  = 11;    tIndices2d3d(0,89)  = 15;// Face 57
    tIndices2d3d(0,90)  = 12;    tIndices2d3d(0,91)  = 19;// Face 58
    tIndices2d3d(0,92)  = 14;    tIndices2d3d(0,93)  = 24;// Face 59
    tIndices2d3d(0,94)  = 16;    tIndices2d3d(0,95)  = 20;// Face 60

    //Offsets
    tOffsets2d1d(0,0)  = 0  ;    tOffsets2d1d(0,1)  = 1  ;    tOffsets2d1d(0,2)  = 2  ;    tOffsets2d1d(0,3)  = 3  ;
    tOffsets2d1d(0,4)  = 4  ;    tOffsets2d1d(0,5)  = 5  ;    tOffsets2d1d(0,6)  = 6  ;    tOffsets2d1d(0,7)  = 7  ;
    tOffsets2d1d(0,8)  = 8  ;    tOffsets2d1d(0,9)  = 9  ;    tOffsets2d1d(0,10) = 10 ;    tOffsets2d1d(0,11) = 11 ;
    tOffsets2d1d(0,12) = 12 ;    tOffsets2d1d(0,13) = 13 ;    tOffsets2d1d(0,14) = 14 ;    tOffsets2d1d(0,15) = 15 ;
    tOffsets2d1d(0,16) = 16 ;    tOffsets2d1d(0,17) = 17 ;    tOffsets2d1d(0,18) = 18 ;    tOffsets2d1d(0,19) = 19 ;
    tOffsets2d1d(0,20) = 20 ;    tOffsets2d1d(0,21) = 21 ;    tOffsets2d1d(0,22) = 22 ;    tOffsets2d1d(0,23) = 23 ;
    tOffsets2d1d(0,24) = 24 ;    tOffsets2d1d(0,25) = 26 ;    tOffsets2d1d(0,26) = 28 ;    tOffsets2d1d(0,27) = 30 ;
    tOffsets2d1d(0,28) = 32 ;    tOffsets2d1d(0,29) = 34 ;    tOffsets2d1d(0,30) = 36 ;    tOffsets2d1d(0,31) = 38 ;
    tOffsets2d1d(0,32) = 40 ;    tOffsets2d1d(0,33) = 42 ;    tOffsets2d1d(0,34) = 44 ;    tOffsets2d1d(0,35) = 46 ;
    tOffsets2d1d(0,36) = 48 ;    tOffsets2d1d(0,37) = 50 ;    tOffsets2d1d(0,38) = 52 ;    tOffsets2d1d(0,39) = 54 ;
    tOffsets2d1d(0,40) = 56 ;    tOffsets2d1d(0,41) = 58 ;    tOffsets2d1d(0,42) = 60 ;    tOffsets2d1d(0,43) = 62 ;
    tOffsets2d1d(0,44) = 64 ;    tOffsets2d1d(0,45) = 66 ;    tOffsets2d1d(0,46) = 68 ;    tOffsets2d1d(0,47) = 70 ;
    tOffsets2d1d(0,48) = 72 ;    tOffsets2d1d(0,49) = 74 ;    tOffsets2d1d(0,50) = 76 ;    tOffsets2d1d(0,51) = 78 ;
    tOffsets2d1d(0,52) = 80 ;    tOffsets2d1d(0,53) = 82 ;    tOffsets2d1d(0,54) = 84 ;    tOffsets2d1d(0,55) = 86 ;
    tOffsets2d1d(0,56) = 88 ;    tOffsets2d1d(0,57) = 90 ;    tOffsets2d1d(0,58) = 92 ;    tOffsets2d1d(0,59) = 94 ;
    tOffsets2d1d(0,60) = 96 ;