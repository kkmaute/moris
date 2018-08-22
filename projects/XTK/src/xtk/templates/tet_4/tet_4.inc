
std::shared_ptr<Matrix<Integer>> tIndices0d1d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 2, 3, 1,  4, 0, 5, 1, 2, 3, 4, 5}});
std::shared_ptr<Matrix<Integer>> tIndices0d2d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 2, 3, 1, 0, 3, 1, 2, 3, 0, 1, 2}});
std::shared_ptr<Matrix<Integer>> tIndices0d3d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 0, 0, 0}});
std::shared_ptr<Matrix<Integer>> tIndices1d0d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 1,  1, 2, 0, 2, 0, 3,  1, 3,  2, 3}});
std::shared_ptr<Matrix<Integer>> tIndices1d2d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 3, 1, 3, 2, 3, 0, 2, 0, 1, 1, 2}});
std::shared_ptr<Matrix<Integer>> tIndices1d3d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 0, 0, 0, 0, 0}});
std::shared_ptr<Matrix<Integer>> tIndices2d0d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 1, 3,
    1, 2, 3,
    0, 2, 3,
    0, 1, 2}});
std::shared_ptr<Matrix<Integer>> tIndices2d1d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 4, 3, 1, 5, 4, 2, 5, 3, 0, 1, 2}});
std::shared_ptr<Matrix<Integer>> tIndices2d3d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 0, 0, 0}});
std::shared_ptr<Matrix<Integer>> tIndices3d0d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 1, 2, 3}});
std::shared_ptr<Matrix<Integer>> tIndices3d1d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 1, 2, 3, 4, 5}});
std::shared_ptr<Matrix<Integer>> tIndices3d2d = mMatrixFactory.create_integer_type_matrix_base({{ 0, 1, 2, 3}});
std::shared_ptr<Matrix<Integer>> tOffsets0d1d = mMatrixFactory.create_integer_type_matrix_base({{0, 3, 6, 9, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets0d2d = mMatrixFactory.create_integer_type_matrix_base({{0, 3, 6, 9, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets0d3d = mMatrixFactory.create_integer_type_matrix_base({{0, 1, 2, 3, 4}});
std::shared_ptr<Matrix<Integer>> tOffsets1d0d = mMatrixFactory.create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets1d2d = mMatrixFactory.create_integer_type_matrix_base({{0, 2, 4, 6, 8, 10, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets1d3d = mMatrixFactory.create_integer_type_matrix_base({{0, 1, 2, 3, 4, 5, 6}});
std::shared_ptr<Matrix<Integer>> tOffsets2d0d = mMatrixFactory.create_integer_type_matrix_base({{0, 3, 6, 9, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets2d1d = mMatrixFactory.create_integer_type_matrix_base({{0, 3, 6, 9, 12}});
std::shared_ptr<Matrix<Integer>> tOffsets2d3d = mMatrixFactory.create_integer_type_matrix_base({{0, 1, 2, 3, 4}});
std::shared_ptr<Matrix<Integer>> tOffsets3d0d = mMatrixFactory.create_integer_type_matrix_base({{0, 4}});
std::shared_ptr<Matrix<Integer>> tOffsets3d1d = mMatrixFactory.create_integer_type_matrix_base({{0, 6}});
std::shared_ptr<Matrix<Integer>> tOffsets3d2d = mMatrixFactory.create_integer_type_matrix_base({{0, 4}});

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceIndices1D    = mMatrixFactory->create_integer_type_matrix_base({{0,1,2,3,4,5}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceRanks1D  = mMatrixFactory->create_integer_type_matrix_base({{1,1,1,1,1,1}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceIndices2D    = mMatrixFactory->create_integer_type_matrix_base({{0,1,2,3}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceRanks2D  = mMatrixFactory->create_integer_type_matrix_base({{2,2,2,2}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceIndices3D  = mMatrixFactory->create_integer_type_matrix_base({{0}});
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tInheritanceRanks3D  = mMatrixFactory->create_integer_type_matrix_base({{3}});




