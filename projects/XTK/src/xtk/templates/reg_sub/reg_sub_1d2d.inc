    // 1D to 2D ---------------------------------------------------------------------------------------------------------------
    std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d2d = mMatrixFactory->create_integer_type_matrix_base(1,180);
    std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d2d = mMatrixFactory->create_integer_type_matrix_base(1,51);

    // Edge 1
    (*tIndices1d2d)(0,0)   =  10 ;    (*tIndices1d2d)(0,1)   =  23 ;    (*tIndices1d2d)(0,2)   =  56 ;

    // Edge 2
    (*tIndices1d2d)(0,3)   =  6  ;    (*tIndices1d2d)(0,4)   =  22 ;    (*tIndices1d2d)(0,5)   =  53 ;

    // Edge 3
    (*tIndices1d2d)(0,6)   =  2  ;    (*tIndices1d2d)(0,7)   =  21 ;    (*tIndices1d2d)(0,8)   =  50 ;

    // Edge 4
    (*tIndices1d2d)(0,9)   =  14 ;    (*tIndices1d2d)(0,10)  =  24 ;    (*tIndices1d2d)(0,11)  =  59 ;

    // Edge 5
    (*tIndices1d2d)(0,12)  =  12 ;    (*tIndices1d2d)(0,13)  =  19 ;    (*tIndices1d2d)(0,14)  =  58 ;

    // Edge 6
    (*tIndices1d2d)(0,15)  =  8  ;    (*tIndices1d2d)(0,16)  =  18 ;    (*tIndices1d2d)(0,17)  =  55 ;

    // Edge 7
    (*tIndices1d2d)(0,18)  =  4  ;    (*tIndices1d2d)(0,19)  =  17 ;    (*tIndices1d2d)(0,20)  =  52 ;

    // Edge 8
    (*tIndices1d2d)(0,21)  =  16 ;    (*tIndices1d2d)(0,22)  =  20 ;    (*tIndices1d2d)(0,23)  =  60 ;

    // Edge 9
    (*tIndices1d2d)(0,24)  =  11 ;    (*tIndices1d2d)(0,25)  =  15 ;    (*tIndices1d2d)(0,26)  =  57 ;

    // Edge 10
    (*tIndices1d2d)(0,27)  =  7  ;    (*tIndices1d2d)(0,28)  =  9  ;    (*tIndices1d2d)(0,29)  =  54 ;

    // Edge 11
    (*tIndices1d2d)(0,30)  =  1  ;    (*tIndices1d2d)(0,31)  =  5  ;    (*tIndices1d2d)(0,32)  =  49 ;

    // Edge 12
    (*tIndices1d2d)(0,33)  =  3  ;    (*tIndices1d2d)(0,34)  =  13 ;    (*tIndices1d2d)(0,35)  =  51 ;

    // Edge 13
    (*tIndices1d2d)(0,36)  =  26 ;    (*tIndices1d2d)(0,37)  =  37 ;    (*tIndices1d2d)(0,38)  =  48 ;
    (*tIndices1d2d)(0,39)  =  50 ;    (*tIndices1d2d)(0,40)  =  51 ;    (*tIndices1d2d)(0,41)  =  59 ;

    // Edge 14
    (*tIndices1d2d)(0,42)  =  27 ;    (*tIndices1d2d)(0,43)  =  40 ;    (*tIndices1d2d)(0,44)  =  44 ;
    (*tIndices1d2d)(0,45)  =  51 ;    (*tIndices1d2d)(0,46)  =  52 ;    (*tIndices1d2d)(0,47)  =  60 ;

    // Edge 15
    (*tIndices1d2d)(0,48)  =  28 ;    (*tIndices1d2d)(0,49)  =  32 ;    (*tIndices1d2d)(0,50)  =  41 ;
    (*tIndices1d2d)(0,51)  =  49 ;    (*tIndices1d2d)(0,52)  =  52 ;    (*tIndices1d2d)(0,53)  =  55 ;

    // Edge 16
    (*tIndices1d2d)(0,54)  =  25 ;    (*tIndices1d2d)(0,55)  =  29 ;    (*tIndices1d2d)(0,56)  =  45 ;
    (*tIndices1d2d)(0,57)  =  49 ;    (*tIndices1d2d)(0,58)  =  50 ;    (*tIndices1d2d)(0,59)  =  53 ;

    // Edge 17
    (*tIndices1d2d)(0,60)  =  31 ;    (*tIndices1d2d)(0,61)  =  36 ;    (*tIndices1d2d)(0,62)  =  42 ;
    (*tIndices1d2d)(0,63)  =  54 ;    (*tIndices1d2d)(0,64)  =  55 ;    (*tIndices1d2d)(0,65)  =  58 ;

    // Edge 18
    (*tIndices1d2d)(0,66)  =  30 ;    (*tIndices1d2d)(0,67)  =  33 ;    (*tIndices1d2d)(0,68)  =  46 ;
    (*tIndices1d2d)(0,69)  =  53 ;    (*tIndices1d2d)(0,70)  =  54 ;    (*tIndices1d2d)(0,71)  =  56 ;

    // Edge 19
    (*tIndices1d2d)(0,72)  =  34 ;    (*tIndices1d2d)(0,73)  =  38 ;    (*tIndices1d2d)(0,74)  =  47 ;
    (*tIndices1d2d)(0,75)  =  56 ;    (*tIndices1d2d)(0,76)  =  57 ;    (*tIndices1d2d)(0,77)  =  59 ;

    // Edge 20
    (*tIndices1d2d)(0,78)  =  35 ;    (*tIndices1d2d)(0,79)  =  39 ;    (*tIndices1d2d)(0,80)  =  43 ;
    (*tIndices1d2d)(0,81)  =  57 ;    (*tIndices1d2d)(0,82)  =  58 ;    (*tIndices1d2d)(0,83)  =  60 ;

    // Edge 21
    (*tIndices1d2d)(0,84)  =  2  ;    (*tIndices1d2d)(0,85)  =  3  ;    (*tIndices1d2d)(0,86)  =  26 ;

    // Edge 22
    (*tIndices1d2d)(0,87)  =  3  ;    (*tIndices1d2d)(0,88)  =  4  ;    (*tIndices1d2d)(0,89)  =  27 ;

    // Edge 23
    (*tIndices1d2d)(0,90)  =  1  ;    (*tIndices1d2d)(0,91)  =  4  ;    (*tIndices1d2d)(0,92)  =  28 ;

    // Edge 24
    (*tIndices1d2d)(0,93)  =  1  ;    (*tIndices1d2d)(0,94)  =  2  ;    (*tIndices1d2d)(0,95)  =  25 ;

    // Edge 25
    (*tIndices1d2d)(0,96)  =  5  ;    (*tIndices1d2d)(0,97)  =  8  ;    (*tIndices1d2d)(0,98)  =  32 ;

    // Edge 26
    (*tIndices1d2d)(0,99)  =  5  ;    (*tIndices1d2d)(0,100) =  6  ;    (*tIndices1d2d)(0,101) =  29 ;

    // Edge 27
    (*tIndices1d2d)(0,102) =  6  ;    (*tIndices1d2d)(0,103) =  7  ;    (*tIndices1d2d)(0,104) =  30 ;

    // Edge 28
    (*tIndices1d2d)(0,105) =  7  ;    (*tIndices1d2d)(0,106) =  8  ;    (*tIndices1d2d)(0,107) =  31 ;

    // Edge 29
    (*tIndices1d2d)(0,108) =  11 ;    (*tIndices1d2d)(0,109) =  12 ;    (*tIndices1d2d)(0,110) =  35 ;

    // Edge 30
    (*tIndices1d2d)(0,111) =  9  ;    (*tIndices1d2d)(0,112) =  12 ;    (*tIndices1d2d)(0,113) =  36 ;

    // Edge 31
    (*tIndices1d2d)(0,114) =  9  ;    (*tIndices1d2d)(0,115) =  10 ;    (*tIndices1d2d)(0,116) =  33 ;

    // Edge 32
    (*tIndices1d2d)(0,117) =  10 ;    (*tIndices1d2d)(0,118) =  11 ;    (*tIndices1d2d)(0,119) =  34 ;

    // Edge 33
    (*tIndices1d2d)(0,120) =  13 ;    (*tIndices1d2d)(0,121) =  16 ;    (*tIndices1d2d)(0,122) =  40 ;

    // Edge 34
    (*tIndices1d2d)(0,123) =  13 ;    (*tIndices1d2d)(0,124) =  14 ;    (*tIndices1d2d)(0,125) =  37 ;

    // Edge 35
    (*tIndices1d2d)(0,126) =  14 ;    (*tIndices1d2d)(0,127) =  15 ;    (*tIndices1d2d)(0,128) =  38 ;

    // Edge 36
    (*tIndices1d2d)(0,129) =  15 ;    (*tIndices1d2d)(0,130) =  16 ;    (*tIndices1d2d)(0,131) =  39 ;

    // Edge 37
    (*tIndices1d2d)(0,132) =  17 ;    (*tIndices1d2d)(0,133) =  20 ;    (*tIndices1d2d)(0,134) =  44 ;

    // Edge 38
    (*tIndices1d2d)(0,135) =  17 ;    (*tIndices1d2d)(0,136) =  18 ;    (*tIndices1d2d)(0,137) =  41 ;

    // Edge 39
    (*tIndices1d2d)(0,138) =  18 ;    (*tIndices1d2d)(0,139) =  19 ;    (*tIndices1d2d)(0,140) =  42 ;

    // Edge 40
    (*tIndices1d2d)(0,141) =  19 ;    (*tIndices1d2d)(0,142) =  20 ;    (*tIndices1d2d)(0,143) =  43 ;

    // Edge 41
    (*tIndices1d2d)(0,144) =  21 ;    (*tIndices1d2d)(0,145) =  24 ;    (*tIndices1d2d)(0,146) =  48 ;

    // Edge 42
    (*tIndices1d2d)(0,147) =  21 ;    (*tIndices1d2d)(0,148) =  22 ;    (*tIndices1d2d)(0,149) =  45 ;

    // Edge 43
    (*tIndices1d2d)(0,150) =  22 ;    (*tIndices1d2d)(0,151) =  23 ;    (*tIndices1d2d)(0,152) =  46 ;

    // Edge 44
    (*tIndices1d2d)(0,153) =  23 ;    (*tIndices1d2d)(0,154) =  24 ;    (*tIndices1d2d)(0,155) =  47 ;

    // Edge 45
    (*tIndices1d2d)(0,156) =  25 ;    (*tIndices1d2d)(0,157) =  26 ;    (*tIndices1d2d)(0,158) =  27 ;
    (*tIndices1d2d)(0,159) =  28 ;

    // Edge 46
    (*tIndices1d2d)(0,160) =  29 ;    (*tIndices1d2d)(0,161) =  30 ;    (*tIndices1d2d)(0,162) =  31 ;
    (*tIndices1d2d)(0,163) =  32 ;

    // Edge 47
    (*tIndices1d2d)(0,164) =  33 ;    (*tIndices1d2d)(0,165) =  34 ;    (*tIndices1d2d)(0,166) =  35 ;
    (*tIndices1d2d)(0,167) =  36 ;

    // Edge 48
    (*tIndices1d2d)(0,168) =  37 ;    (*tIndices1d2d)(0,169) =  38 ;    (*tIndices1d2d)(0,170) =  39 ;
    (*tIndices1d2d)(0,171) =  40 ;

    // Edge 49
    (*tIndices1d2d)(0,172) =  41 ;    (*tIndices1d2d)(0,173) =  42 ;    (*tIndices1d2d)(0,174) =  43 ;
    (*tIndices1d2d)(0,175) =  44 ;

    // Edge 50
    (*tIndices1d2d)(0,176) =  45 ;    (*tIndices1d2d)(0,177) =  46 ;    (*tIndices1d2d)(0,178) =  47 ;
    (*tIndices1d2d)(0,179) =  48 ;

    //Offsets
    (*tOffsets1d2d)(0,0)   = 0  ;    (*tOffsets1d2d)(0,1)   = 3  ;    (*tOffsets1d2d)(0,2)   = 6  ;    (*tOffsets1d2d)(0,3)   = 9;
    (*tOffsets1d2d)(0,4)   = 12 ;    (*tOffsets1d2d)(0,5)   = 15 ;    (*tOffsets1d2d)(0,6)   = 18 ;    (*tOffsets1d2d)(0,7)   = 21;
    (*tOffsets1d2d)(0,8)   = 24 ;    (*tOffsets1d2d)(0,9)   = 27 ;    (*tOffsets1d2d)(0,10)  = 30 ;    (*tOffsets1d2d)(0,11)  = 33;
    (*tOffsets1d2d)(0,12)  = 36 ;    (*tOffsets1d2d)(0,13)  = 42 ;    (*tOffsets1d2d)(0,14)  = 48 ;    (*tOffsets1d2d)(0,15)  = 54;
    (*tOffsets1d2d)(0,16)  = 60 ;    (*tOffsets1d2d)(0,17)  = 66 ;    (*tOffsets1d2d)(0,18)  = 72 ;    (*tOffsets1d2d)(0,19)  = 78;
    (*tOffsets1d2d)(0,20)  = 84 ;    (*tOffsets1d2d)(0,21)  = 87 ;    (*tOffsets1d2d)(0,22)  = 90 ;    (*tOffsets1d2d)(0,23)  = 93;
    (*tOffsets1d2d)(0,24)  = 96 ;    (*tOffsets1d2d)(0,25)  = 99 ;    (*tOffsets1d2d)(0,26)  = 102;    (*tOffsets1d2d)(0,27)  = 105;
    (*tOffsets1d2d)(0,28)  = 108;    (*tOffsets1d2d)(0,29)  = 111;    (*tOffsets1d2d)(0,30)  = 114;    (*tOffsets1d2d)(0,31)  = 117;
    (*tOffsets1d2d)(0,32)  = 120;    (*tOffsets1d2d)(0,33)  = 123;    (*tOffsets1d2d)(0,34)  = 126;    (*tOffsets1d2d)(0,35)  = 129;
    (*tOffsets1d2d)(0,36)  = 132;    (*tOffsets1d2d)(0,37)  = 135;    (*tOffsets1d2d)(0,38)  = 138;    (*tOffsets1d2d)(0,39)  = 141;
    (*tOffsets1d2d)(0,40)  = 144;    (*tOffsets1d2d)(0,41)  = 147;    (*tOffsets1d2d)(0,42)  = 150;    (*tOffsets1d2d)(0,43)  = 153;
    (*tOffsets1d2d)(0,44)  = 156;    (*tOffsets1d2d)(0,45)  = 160;    (*tOffsets1d2d)(0,46)  = 164;    (*tOffsets1d2d)(0,47)  = 168;
    (*tOffsets1d2d)(0,48)  = 172;    (*tOffsets1d2d)(0,49)  = 176;    (*tOffsets1d2d)(0,50)  = 180;

