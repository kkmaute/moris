    // 0D to 3D ---------------------------------------------------------------------------------------------------------------
    std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d3d = mMatrixFactory->create_integer_type_matrix_base(1,96);
    std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d3d = mMatrixFactory->create_integer_type_matrix_base(1,16);

    //Node 1 (n1)
    (*tIndices0d3d)(0,0)   = 10;    (*tIndices0d3d)(0,1)   = 11;    (*tIndices0d3d)(0,2)   = 14;    (*tIndices0d3d)(0,3)   = 15;
    (*tIndices0d3d)(0,4)   = 23;    (*tIndices0d3d)(0,5)   = 24;

    //Node 2 (n2)
    (*tIndices0d3d)(0,6)   = 6 ;    (*tIndices0d3d)(0,7)   = 7 ;    (*tIndices0d3d)(0,8)   = 9 ;    (*tIndices0d3d)(0,9)   = 10;
    (*tIndices0d3d)(0,10)  = 22;    (*tIndices0d3d)(0,11)  = 23;

    //Node 3 (n3)
    (*tIndices0d3d)(0,12)  = 1 ;    (*tIndices0d3d)(0,13)  = 2 ;    (*tIndices0d3d)(0,14)  = 5 ;    (*tIndices0d3d)(0,15)  = 6;
    (*tIndices0d3d)(0,16)  = 21;    (*tIndices0d3d)(0,17)  = 22;

    //Node 4 (n4)
    (*tIndices0d3d)(0,18)  = 2 ;    (*tIndices0d3d)(0,19)  = 3 ;    (*tIndices0d3d)(0,20)  = 13;    (*tIndices0d3d)(0,21)  = 14;
    (*tIndices0d3d)(0,22)  = 21;    (*tIndices0d3d)(0,23)  = 24;

    //Node 5 (n5)
    (*tIndices0d3d)(0,24)  = 11;    (*tIndices0d3d)(0,25)  = 12;    (*tIndices0d3d)(0,26)  = 15;    (*tIndices0d3d)(0,27)  = 16;
    (*tIndices0d3d)(0,28)  = 19;    (*tIndices0d3d)(0,29)  = 20;

    //Node 6 (n6)
    (*tIndices0d3d)(0,30)  = 7 ;    (*tIndices0d3d)(0,31)  = 8 ;    (*tIndices0d3d)(0,32)  = 9 ;    (*tIndices0d3d)(0,33)  = 12;
    (*tIndices0d3d)(0,34)  = 18;    (*tIndices0d3d)(0,35)  = 19;

    //Node 7 (n7)
    (*tIndices0d3d)(0,36)  = 1 ;    (*tIndices0d3d)(0,37)  = 4 ;    (*tIndices0d3d)(0,38)  = 5 ;    (*tIndices0d3d)(0,39)  = 8;
    (*tIndices0d3d)(0,40)  = 17;    (*tIndices0d3d)(0,41)  = 18;

    //Node 8 (n8)
    (*tIndices0d3d)(0,42)  = 3 ;    (*tIndices0d3d)(0,43)  = 4 ;    (*tIndices0d3d)(0,44)  = 13 ;    (*tIndices0d3d)(0,45)  = 16;
    (*tIndices0d3d)(0,46)  = 17;    (*tIndices0d3d)(0,47)  = 20;

    //Node 9 (n9)
    (*tIndices0d3d)(0,48)  = 9 ;    (*tIndices0d3d)(0,49)  = 10;    (*tIndices0d3d)(0,50)  = 11;    (*tIndices0d3d)(0,51)  = 12;

    //Node 10 (n10)
    (*tIndices0d3d)(0,52)  = 5 ;    (*tIndices0d3d)(0,53)  = 6 ;    (*tIndices0d3d)(0,54)  = 7 ;    (*tIndices0d3d)(0,55)  = 8;

    //Node 11 (n11)
    (*tIndices0d3d)(0,56)  = 1 ;    (*tIndices0d3d)(0,57)  = 2 ;    (*tIndices0d3d)(0,58)  = 3 ;    (*tIndices0d3d)(0,59)  = 4;

    //Node 12 (n12)
    (*tIndices0d3d)(0,60)  = 13;    (*tIndices0d3d)(0,61)  = 14;    (*tIndices0d3d)(0,62)  = 15;    (*tIndices0d3d)(0,63)  = 16;

    //Node 13 (n13)
    (*tIndices0d3d)(0,64)  = 21;    (*tIndices0d3d)(0,65)  = 22;    (*tIndices0d3d)(0,66)  = 23;    (*tIndices0d3d)(0,67)  = 24;

    //Node 14 (n14)
    (*tIndices0d3d)(0,68)  = 17;    (*tIndices0d3d)(0,69)  = 18;    (*tIndices0d3d)(0,70)  = 19;    (*tIndices0d3d)(0,71)  = 20;

    //Node 15 (n15)
    (*tIndices0d3d)(0,72)  = 1 ;    (*tIndices0d3d)(0,73)  = 2 ;    (*tIndices0d3d)(0,74)  = 3 ;    (*tIndices0d3d)(0,75)  = 4 ;
    (*tIndices0d3d)(0,76)  = 5 ;    (*tIndices0d3d)(0,77)  = 6 ;    (*tIndices0d3d)(0,78)  = 7 ;    (*tIndices0d3d)(0,79)  = 8 ;
    (*tIndices0d3d)(0,80)  = 9 ;    (*tIndices0d3d)(0,81)  = 10;    (*tIndices0d3d)(0,82)  = 11;    (*tIndices0d3d)(0,83)  = 12;
    (*tIndices0d3d)(0,84)  = 13;    (*tIndices0d3d)(0,85)  = 14;    (*tIndices0d3d)(0,86)  = 15;    (*tIndices0d3d)(0,87)  = 16;
    (*tIndices0d3d)(0,88)  = 17;    (*tIndices0d3d)(0,89)  = 18;    (*tIndices0d3d)(0,90)  = 19;    (*tIndices0d3d)(0,91)  = 20;
    (*tIndices0d3d)(0,92)  = 21;    (*tIndices0d3d)(0,93)  = 22;    (*tIndices0d3d)(0,94)  = 23;    (*tIndices0d3d)(0,95)  = 24;


    //Offsets
    (*tOffsets0d3d)(0,0)   = 0 ;    (*tOffsets0d3d)(0,1)   = 6 ;    (*tOffsets0d3d)(0,2)   = 12;    (*tOffsets0d3d)(0,3)   = 18;
    (*tOffsets0d3d)(0,4)   = 24;    (*tOffsets0d3d)(0,5)   = 30;    (*tOffsets0d3d)(0,6)   = 36;    (*tOffsets0d3d)(0,7)   = 42;
    (*tOffsets0d3d)(0,8)   = 48;    (*tOffsets0d3d)(0,9)   = 52;    (*tOffsets0d3d)(0,10)  = 56;    (*tOffsets0d3d)(0,11)  = 60;
    (*tOffsets0d3d)(0,12)  = 64;    (*tOffsets0d3d)(0,13)  = 68;    (*tOffsets0d3d)(0,14)  = 72;    (*tOffsets0d3d)(0,15)  = 96;
