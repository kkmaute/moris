
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 6, 12, 13, 14, 15, 16, 17, 18, 15, 7, 8, 9, 6, 1, 4, 0, 18, 10, 11, 9, 12, 5, 1, 2, 7, 10, 13, 16, 3, 4, 5, 8, 11, 14, 17}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 4, 6, 10, 11, 12, 13, 15, 14, 16, 17, 13, 15, 5, 7, 9, 4, 6, 1, 0, 3, 16, 17, 8, 7, 9, 10, 12, 1, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15, 17, 0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 3, 4, 5, 4, 2, 1, 0, 5, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 4, 6, 0, 6, 0, 7, 4, 7, 6, 7, 0, 3, 3, 6, 3, 7, 3, 5, 5, 6, 5, 7, 1, 5, 1, 6, 1, 7, 1, 2, 2, 6, 2, 7, 2, 4}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 1, 3, 17, 2, 3, 6, 0, 2, 4, 0, 16, 1, 1, 2, 5, 8, 11, 14, 4, 6, 5, 6, 9, 4, 5, 7, 7, 9, 8, 9, 12, 7, 10, 8, 10, 12, 11, 12, 15, 10, 11, 13, 13, 15, 14, 15, 17, 13, 14, 16, 16, 17}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 5, 0, 1, 0, 1, 0, 5, 0, 1, 2, 3, 4, 5, 1, 1, 2, 2, 1, 2, 2, 3, 2, 3, 3, 3, 4, 3, 4, 4, 4, 5, 5, 4, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 7, 4, 6, 7, 0, 6, 7, 0, 4, 6, 0, 3, 7, 3, 6, 7, 0, 3, 6, 3, 5, 7, 5, 6, 7, 3, 5, 6, 1, 5, 7, 1, 6, 7, 1, 5, 6, 1, 2, 7, 2, 6, 7, 1, 2, 6, 2, 4, 7, 2, 4, 6}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 3, 1, 5, 4, 2, 5, 3, 0, 1, 2, 6, 8, 3, 7, 5, 8, 6, 7, 2, 9, 11, 8, 10, 5, 11, 9, 10, 7, 12, 11, 14, 13, 5, 14, 12, 10, 13, 15, 17, 14, 16, 5, 17, 15, 16, 13, 18, 4, 17, 18, 1, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 5, 0, 1, 0, 1, 1, 2, 1, 2, 2, 3, 2, 3, 3, 4, 3, 4, 4, 5, 4, 5, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d0d = mMatrixFactory->create_integer_type_matrix_base({{ 4, 0, 6, 7,
    0, 3, 6, 7,
    3, 5, 6, 7,
    5, 1, 6, 7,
    1, 2, 6, 7,
    2, 4, 6, 7}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5, 6, 7, 2, 3, 8, 5, 9, 10, 7, 8, 11, 5, 12, 10, 13, 14, 11, 5, 15, 16, 13, 14, 17, 5, 18, 1, 16, 17, 4, 5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5, 2, 6, 7, 8, 5, 9, 10, 8, 11, 12, 13, 14, 11, 15, 16, 1, 14, 17}});

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16, 20, 24, 31, 38}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 5, 10, 15, 20, 25, 30, 42, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12, 18, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 5, 8, 11, 14, 20, 22, 25, 28, 30, 33, 36, 38, 41, 44, 46, 49, 52, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 3, 5, 7, 9, 15, 16, 18, 20, 21, 23, 25, 26, 28, 30, 31, 33, 35, 36}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16, 20, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 6, 12, 18, 24, 30, 36}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16, 20, 24}});
