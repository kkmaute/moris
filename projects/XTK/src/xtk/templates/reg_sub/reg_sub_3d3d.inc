    // 3D to 3D ---------------------------------------------------------------------------------------------------------------
    // Required elements share a face to be considered neighbors
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d3d = mMatrixFactory->create_integer_type_matrix_base(1, 144);
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d3d = mMatrixFactory->create_integer_type_matrix_base(1, 25);

    // Subelement 1
    (*tIndices3d3d)(0,0)   =  2 ;    (*tIndices3d3d)(0,1)   =  4 ;    (*tIndices3d3d)(0,2)   =  5 ;

    // Subelement 2
    (*tIndices3d3d)(0,3)   =  1 ;    (*tIndices3d3d)(0,4)   =  3 ;    (*tIndices3d3d)(0,5)   =  21;

    // Subelement 3
    (*tIndices3d3d)(0,6)   =  2 ;    (*tIndices3d3d)(0,7)   =  4 ;    (*tIndices3d3d)(0,8)   =  13;

    // Subelement 4
    (*tIndices3d3d)(0,9)   =  1 ;    (*tIndices3d3d)(0,10)  =  3 ;    (*tIndices3d3d)(0,11)  =  17;

    // Subelement 5
    (*tIndices3d3d)(0,12)  =  1 ;    (*tIndices3d3d)(0,13)  =  6 ;    (*tIndices3d3d)(0,14)  =  8 ;

    // Subelement 6
    (*tIndices3d3d)(0,15)  =  5 ;    (*tIndices3d3d)(0,16)  =  7 ;    (*tIndices3d3d)(0,17)  =  22;

    // Subelement 7
    (*tIndices3d3d)(0,18)  =  6 ;    (*tIndices3d3d)(0,19)  =  8 ;    (*tIndices3d3d)(0,20)  =  9 ;

    // Subelement 8
    (*tIndices3d3d)(0,21)  =  5 ;    (*tIndices3d3d)(0,22)  =  7 ;    (*tIndices3d3d)(0,23)  =  18;

    // Subelement 9
    (*tIndices3d3d)(0,24)  =  7 ;    (*tIndices3d3d)(0,25)  =  10;    (*tIndices3d3d)(0,26)  =  12;

    // Subelement 10
    (*tIndices3d3d)(0,27)  =  9 ;    (*tIndices3d3d)(0,28)  =  11;    (*tIndices3d3d)(0,29)  =  23;

    // Subelement 11
    (*tIndices3d3d)(0,30)  =  10;    (*tIndices3d3d)(0,31)  =  12;    (*tIndices3d3d)(0,32)  =  15;

    // Subelement 12
    (*tIndices3d3d)(0,33)  =  9 ;    (*tIndices3d3d)(0,34)  =  11;    (*tIndices3d3d)(0,35)  =  19;

    // Subelement 13
    (*tIndices3d3d)(0,36)  =  3 ;    (*tIndices3d3d)(0,37)  =  14;    (*tIndices3d3d)(0,38)  =  16;

    // Subelement 14
    (*tIndices3d3d)(0,39)  =  13;    (*tIndices3d3d)(0,40)  =  15;    (*tIndices3d3d)(0,41)  =  24;

    // Subelement 15
    (*tIndices3d3d)(0,42)  =  11;    (*tIndices3d3d)(0,43)  =  14;    (*tIndices3d3d)(0,44)  =  16;

    // Subelement 16
    (*tIndices3d3d)(0,45)  =  13;    (*tIndices3d3d)(0,46)  =  15;    (*tIndices3d3d)(0,47)  =  20;

    // Subelement 17
    (*tIndices3d3d)(0,48)  =  4 ;    (*tIndices3d3d)(0,49)  =  18;    (*tIndices3d3d)(0,50)  =  20;

    // Subelement 18
    (*tIndices3d3d)(0,51)  =  8 ;    (*tIndices3d3d)(0,52)  =  17;    (*tIndices3d3d)(0,53)  =  19;

    // Subelement 19
    (*tIndices3d3d)(0,54)  =  12;    (*tIndices3d3d)(0,55)  =  18;    (*tIndices3d3d)(0,56)  =  20;

    // Subelement 20
    (*tIndices3d3d)(0,57)  =  16;    (*tIndices3d3d)(0,58)  =  17;    (*tIndices3d3d)(0,59)  =  19;

    // Subelement 21
    (*tIndices3d3d)(0,60)  =  2 ;    (*tIndices3d3d)(0,61)  =  22;    (*tIndices3d3d)(0,62)  =  24;

    // Subelement 22
    (*tIndices3d3d)(0,63)  =  6 ;    (*tIndices3d3d)(0,64)  =  21;    (*tIndices3d3d)(0,65)  =  23;

    // Subelement 23
    (*tIndices3d3d)(0,66)  =  10;    (*tIndices3d3d)(0,67)  =  22;    (*tIndices3d3d)(0,68)  =  24;

    // Subelement 24
    (*tIndices3d3d)(0,69)  =  14;    (*tIndices3d3d)(0,70)  =  21;    (*tIndices3d3d)(0,71)  =  23;

    //Offsets
    (*tOffsets3d3d)(0,0)  = 0 ;    (*tOffsets3d3d)(0,1)  = 3 ;    (*tOffsets3d3d)(0,2)  = 6 ;    (*tOffsets3d3d)(0,3)  = 9 ;
    (*tOffsets3d3d)(0,4)  = 12;    (*tOffsets3d3d)(0,5)  = 15;    (*tOffsets3d3d)(0,6)  = 18;    (*tOffsets3d3d)(0,7)  = 21;
    (*tOffsets3d3d)(0,8)  = 24;    (*tOffsets3d3d)(0,9)  = 27;    (*tOffsets3d3d)(0,10) = 30;    (*tOffsets3d3d)(0,11) = 33;
    (*tOffsets3d3d)(0,12) = 36;    (*tOffsets3d3d)(0,13) = 39;    (*tOffsets3d3d)(0,14) = 42;    (*tOffsets3d3d)(0,15) = 45;
    (*tOffsets3d3d)(0,16) = 48;    (*tOffsets3d3d)(0,17) = 51;    (*tOffsets3d3d)(0,18) = 54;    (*tOffsets3d3d)(0,19) = 57;
    (*tOffsets3d3d)(0,20) = 60;    (*tOffsets3d3d)(0,21) = 63;    (*tOffsets3d3d)(0,22) = 66;    (*tOffsets3d3d)(0,23) = 69;
    (*tOffsets3d3d)(0,24) = 72;
