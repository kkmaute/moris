//moris::Mat<moris::uint> tIndices3d0d(1,24);
//moris::Mat<moris::uint> tOffsets3d0d(1,7);

// Node structure

// High and low are opposite on interface
// {0, 1, 2, 3, l, ml, mh, h}; where i0 means interface node 0 and i0 is the lowest id interface node
// {0, 1, 2, 3, 4,  5,  6, 7}

//// Nodes connected to element
//// Element 1
//tIndices3d0d(0,0) =  0;   tIndices3d0d(0,1) = 3;   tIndices3d0d(0,2) = 7;   tIndices3d0d(0,3) = 6;
//
//// Element 2
//tIndices3d0d(0,4) =  3;   tIndices3d0d(0,5) = 4;   tIndices3d0d(0,6) = 6;   tIndices3d0d(0,7) = 7;
//
//// Element 3
//tIndices3d0d(0,8) =  3;  tIndices3d0d(0,9) = 4;   tIndices3d0d(0,10) = 7;   tIndices3d0d(0,11) = 5;
//
//// Element 4
//tIndices3d0d(0,12) = 1;  tIndices3d0d(0,13) = 2;   tIndices3d0d(0,14) = 7;   tIndices3d0d(0,15) = 5;
//
//// Element 5
//tIndices3d0d(0,16) = 2 ; tIndices3d0d(0,17) = 4;   tIndices3d0d(0,18) = 7;   tIndices3d0d(0,19) = 5;
//
//// Element 6
//tIndices3d0d(0,20) = 2;  tIndices3d0d(0,21) = 4;   tIndices3d0d(0,22) = 7;   tIndices3d0d(0,23) = 6;
//
//
//tOffsets3d0d(0,0) = 0 ;  tOffsets3d0d(0,1) = 4;  tOffsets3d0d(0,2) = 8;
//tOffsets3d0d(0,3) = 12;  tOffsets3d0d(0,4) = 16; tOffsets3d0d(0,5) = 20;
//tOffsets3d0d(0,6) = 24;

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 12, 14, 15, 13, 16, 17, 18, 12, 1, 4, 6, 9, 0, 7, 8, 10, 6, 17, 11, 9, 10, 15, 16, 5, 3, 4, 7, 18, 1, 2, 5, 8, 11, 13, 14}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 10, 12, 13, 11, 14, 15, 16, 17, 10, 13, 1, 4, 6, 7, 9, 0, 3, 5, 8, 4, 6, 7, 14, 15, 16, 8, 9, 11, 12, 7, 10, 14, 1, 2, 5, 17, 0, 6, 16, 1, 2, 3, 4, 5, 8, 9, 11, 12, 13, 15, 17}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 4, 5, 3, 1, 2, 0, 1, 2, 4, 5, 2, 3, 4, 1, 0, 5, 0, 2, 3, 4, 5, 1}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 3, 7, 0, 7, 0, 6, 3, 6, 6, 7, 3, 4, 4, 6, 4, 7, 3, 5, 4, 5, 5, 7, 1, 2, 2, 7, 1, 7, 1, 5, 2, 5, 2, 4, 2, 6}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 3, 1, 4, 9, 2, 3, 2, 0, 1, 0, 6, 1, 2, 5, 17, 4, 6, 7, 5, 6, 16, 4, 15, 5, 8, 9, 7, 8, 7, 14, 8, 9, 11, 12, 10, 13, 13, 11, 15, 17, 12, 13, 12, 10, 11, 10, 14, 14, 15, 16, 17, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 2, 1, 0, 0, 1, 0, 0, 1, 5, 1, 2, 1, 5, 2, 4, 5, 1, 2, 2, 4, 2, 3, 4, 3, 3, 4, 5, 3, 3, 4, 3, 4, 5, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 6, 3, 6, 7, 0, 6, 7, 0, 3, 7, 3, 4, 7, 4, 6, 7, 3, 4, 6, 3, 4, 5, 4, 5, 7, 3, 5, 7, 1, 2, 5, 2, 5, 7, 1, 5, 7, 1, 2, 7, 2, 4, 5, 2, 4, 7, 2, 4, 6, 2, 6, 7}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 3, 4, 5, 1, 3, 5, 2, 0, 1, 2, 6, 8, 1, 7, 5, 8, 6, 7, 4, 6, 10, 9, 10, 11, 8, 9, 11, 1, 12, 16, 15, 16, 11, 13, 15, 11, 14, 12, 13, 14, 17, 10, 16, 17, 8, 13, 17, 7, 18, 18, 5, 13}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 1, 0, 0, 1, 2, 1, 5, 1, 2, 2, 4, 2, 3, 3, 4, 3, 3, 4, 4, 5, 5, 5}});

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 7, 6,
    4, 3, 6, 7,
    4, 3, 7, 5,
    1, 2, 7, 5,
    2, 4, 7, 5,
    4, 2, 7, 6}});

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5, 6, 7, 4, 1, 8, 5, 6, 8, 1, 9, 10, 11, 12, 13, 14, 15, 16, 11, 17, 8, 13, 16, 10, 11, 17, 8, 13, 18, 7, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5, 1, 6, 7, 8, 9, 4, 10, 11, 12, 13, 14, 8, 11, 15, 16, 5, 17, 15}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 11, 16, 21, 26, 31, 38}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 13, 20, 28, 35, 42, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 2, 5, 8, 12, 15, 18, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 6, 8, 10, 13, 17, 20, 23, 27, 29, 32, 36, 38, 42, 44, 46, 49, 52, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 4, 5, 6, 8, 11, 13, 15, 19, 20, 22, 25, 26, 29, 30, 31, 33, 35, 36}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16, 20, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 6, 12, 18, 24, 30, 36}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16, 20, 24}});

