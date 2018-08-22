// 3D to 2D ---------------------------------------------------------------------------------------------------------------

std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d2d = mMatrixFactory->create_integer_type_matrix_base(1, 96);
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d2d = mMatrixFactory->create_integer_type_matrix_base(1, 25);

(*tOffsets3d2d)(0,0) = 0;

// Subelement 1
(*tIndices3d2d)(0,0) = 1; (*tIndices3d2d)(0,1) = 25; (*tIndices3d2d)(0,2) = 28; (*tIndices3d2d)(0,3) = 49;
(*tOffsets3d2d)(0,1) = 4;

// Subelement 2
(*tIndices3d2d)(0,4) = 2; (*tIndices3d2d)(0,5) = 25; (*tIndices3d2d)(0,6) = 26; (*tIndices3d2d)(0,7) = 50;
(*tOffsets3d2d)(0,2) = 8;

// Subelement 3
(*tIndices3d2d)(0,8) = 3; (*tIndices3d2d)(0,9) = 26; (*tIndices3d2d)(0,10) = 27; (*tIndices3d2d)(0,11) = 51;
(*tOffsets3d2d)(0,3) = 12;

// Subelement 4
(*tIndices3d2d)(0,12) = 4; (*tIndices3d2d)(0,13) = 27; (*tIndices3d2d)(0,14) = 28; (*tIndices3d2d)(0,15) = 52;
(*tOffsets3d2d)(0,4) = 16;

// Subelement 5
(*tIndices3d2d)(0,16) = 5; (*tIndices3d2d)(0,17) = 29; (*tIndices3d2d)(0,18) = 32; (*tIndices3d2d)(0,19) = 49;
(*tOffsets3d2d)(0,5) = 20;

// Subelement 6
(*tIndices3d2d)(0,20) = 6; (*tIndices3d2d)(0,21) = 29; (*tIndices3d2d)(0,22) = 30; (*tIndices3d2d)(0,23) = 53;
(*tOffsets3d2d)(0,6) = 24;

// Subelement 7
(*tIndices3d2d)(0,24) = 7; (*tIndices3d2d)(0,25) = 30; (*tIndices3d2d)(0,26) = 31; (*tIndices3d2d)(0,27) = 54;
(*tOffsets3d2d)(0,7) = 28;

// Subelement 8
(*tIndices3d2d)(0,28) = 8; (*tIndices3d2d)(0,29) = 31; (*tIndices3d2d)(0,30) = 32; (*tIndices3d2d)(0,31) = 55;
(*tOffsets3d2d)(0,8) = 32;

// Subelement 9
(*tIndices3d2d)(0,32) = 9; (*tIndices3d2d)(0,33) = 33; (*tIndices3d2d)(0,34) = 36; (*tIndices3d2d)(0,35) = 54;
(*tOffsets3d2d)(0,9) = 36;

// Subelement 10
(*tIndices3d2d)(0,36) = 10; (*tIndices3d2d)(0,37) = 33; (*tIndices3d2d)(0,38) = 34; (*tIndices3d2d)(0,39) = 56;
(*tOffsets3d2d)(0,10) = 40;

// Subelement 11
(*tIndices3d2d)(0,40) = 11; (*tIndices3d2d)(0,41) = 34; (*tIndices3d2d)(0,42) = 35; (*tIndices3d2d)(0,43) = 57;
(*tOffsets3d2d)(0,11) = 44;

// Subelement 12
(*tIndices3d2d)(0,44) = 12; (*tIndices3d2d)(0,45) = 35; (*tIndices3d2d)(0,46) = 36; (*tIndices3d2d)(0,47) = 58;
(*tOffsets3d2d)(0,12) = 48;

// Subelement 13
(*tIndices3d2d)(0,48) = 13; (*tIndices3d2d)(0,49) = 37; (*tIndices3d2d)(0,50) = 40; (*tIndices3d2d)(0,51) = 51;
(*tOffsets3d2d)(0,13) = 52;

// Subelement 14
(*tIndices3d2d)(0,52) = 14; (*tIndices3d2d)(0,53) = 37; (*tIndices3d2d)(0,54) = 38; (*tIndices3d2d)(0,55) = 59;
(*tOffsets3d2d)(0,14) = 56;

// Subelement 15
(*tIndices3d2d)(0,56) = 15; (*tIndices3d2d)(0,57) = 38; (*tIndices3d2d)(0,58) = 39; (*tIndices3d2d)(0,59) = 57;
(*tOffsets3d2d)(0,15) = 60;

// Subelement 16
(*tIndices3d2d)(0,60) = 16; (*tIndices3d2d)(0,61) = 39; (*tIndices3d2d)(0,62) = 40; (*tIndices3d2d)(0,63) = 60;
(*tOffsets3d2d)(0,16) = 64;

// Subelement 17
(*tIndices3d2d)(0,64) = 17; (*tIndices3d2d)(0,65) = 41; (*tIndices3d2d)(0,66) = 44; (*tIndices3d2d)(0,67) = 52;
(*tOffsets3d2d)(0,17) = 68;

// Subelement 18
(*tIndices3d2d)(0,68) = 18; (*tIndices3d2d)(0,69) = 41; (*tIndices3d2d)(0,70) = 42; (*tIndices3d2d)(0,71) = 55;
(*tOffsets3d2d)(0,18) = 72;

// Subelement 19
(*tIndices3d2d)(0,72) = 19; (*tIndices3d2d)(0,73) = 42; (*tIndices3d2d)(0,74) = 43; (*tIndices3d2d)(0,75) = 58;
(*tOffsets3d2d)(0,19) = 76;

// Subelement 20
(*tIndices3d2d)(0,76) = 20; (*tIndices3d2d)(0,77) = 43; (*tIndices3d2d)(0,78) = 44; (*tIndices3d2d)(0,79) = 60;
(*tOffsets3d2d)(0,20) = 80;

// Subelement 21
(*tIndices3d2d)(0,80) = 21; (*tIndices3d2d)(0,81) = 45; (*tIndices3d2d)(0,82) = 48; (*tIndices3d2d)(0,83) = 50;
(*tOffsets3d2d)(0,21) = 84;

// Subelement 22
(*tIndices3d2d)(0,84) = 22; (*tIndices3d2d)(0,85) = 45; (*tIndices3d2d)(0,86) = 46; (*tIndices3d2d)(0,87) = 53;
(*tOffsets3d2d)(0,22) = 88;

// Subelement 23
(*tIndices3d2d)(0,88) = 23; (*tIndices3d2d)(0,89) = 46; (*tIndices3d2d)(0,90) = 47; (*tIndices3d2d)(0,91) = 56;
(*tOffsets3d2d)(0,23) = 92;

// Subelement 24
(*tIndices3d2d)(0,92) = 24; (*tIndices3d2d)(0,93) = 47; (*tIndices3d2d)(0,94) = 48; (*tIndices3d2d)(0,95) = 59;
(*tOffsets3d2d)(0,24) = 96;
