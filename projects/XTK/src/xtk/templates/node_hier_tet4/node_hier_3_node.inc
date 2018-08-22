// Node structure
// {0,1,2,3,l,m,h});
// {0,1,2,3,4, 5 ,6}

// Nodes connected to element
// Element 1 (3 interface nodes and 1 regular node  (0,H,L,M)
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 6, 7, 8, 9, 12, 10, 11, 13, 9, 14, 12, 13, 1, 4, 0, 6, 5, 1, 2, 7, 10, 3, 4, 5, 8, 11, 14}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 2, 3, 4, 5, 6, 7, 9, 10, 12, 8, 11, 7, 9, 12, 10, 11, 12, 1, 0, 3, 4, 6, 1, 2, 5, 8, 3, 6, 9, 0, 1, 2, 4, 5, 7, 8, 10, 11}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices0d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 2, 3, 3, 0, 1, 0, 1, 2, 0, 1, 2, 3}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 4, 5, 0, 5, 0, 6, 4, 6, 5, 6, 1, 4, 1, 5, 1, 6, 1, 2, 2, 5, 2, 6, 1, 3, 2, 3, 3, 6}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 3, 1, 3, 6, 2, 3, 0, 2, 0, 4, 1, 1, 2, 5, 8, 4, 6, 5, 6, 9, 4, 5, 7, 10, 7, 9, 12, 8, 9, 7, 8, 11, 10, 12, 11, 12, 10, 11}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices1d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 1, 0, 0, 0, 1, 0, 1, 2, 1, 1, 2, 1, 2, 3, 2, 3, 2, 2, 3, 3, 3, 3}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 6, 4, 5, 6, 0, 5, 6, 0, 4, 5, 1, 4, 6, 1, 5, 6, 1, 4, 5, 1, 2, 6, 2, 5, 6, 1, 2, 5, 1, 3, 6, 2, 3, 6, 1, 2, 3}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 3, 1, 5, 4, 2, 5, 3, 0, 1, 2, 6, 4, 8, 7, 5, 8, 6, 1, 7, 9, 11, 8, 10, 5, 11, 9, 10, 7, 12, 14, 8, 13, 14, 11, 9, 13, 12}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices2d3d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 0, 1, 0, 0, 1, 1, 2, 1, 2, 3, 2, 2, 3, 3, 3}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d0d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 4, 5, 6,4, 1, 5, 6,1, 2, 5, 6,1, 3, 2, 6}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d1d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5, 6, 1, 7, 8, 4, 5, 9, 10, 7, 8, 11, 5, 12, 13, 9, 8, 14, 11}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d2d = mMatrixFactory->create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 1, 5, 6, 7, 8, 5, 9, 10, 11, 7, 12}});

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 8, 12, 15, 19, 24, 30}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 10, 15, 18, 23, 30, 39}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets0d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 4, 6, 7, 9, 12, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 2, 5, 7, 9, 12, 16, 18, 21, 25, 28, 30, 33, 35, 37, 39}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets1d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 3, 4, 5, 7, 10, 11, 13, 16, 18, 19, 21, 22, 23, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets2d3d = mMatrixFactory->create_integer_type_matrix_base({{0, 1, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d0d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d1d = mMatrixFactory->create_integer_type_matrix_base({{0, 6, 12, 18, 24}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d2d = mMatrixFactory->create_integer_type_matrix_base({{0, 4, 8, 12, 16}});
