// 3D to 1D ---------------------------------------------------------------------------------------------------------------
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tIndices3d1d = mMatrixFactory->create_integer_type_matrix_base(1, 144);
std::shared_ptr<xtk::Matrix_Base<Integer, Integer_Matrix>> tOffsets3d1d = mMatrixFactory->create_integer_type_matrix_base(1, 25);

// Subelement 1
(*tIndices3d1d)(0,0) = 11; (*tIndices3d1d)(0,1) = 15; (*tIndices3d1d)(0,2) = 16;
(*tIndices3d1d)(0,3) = 23; (*tIndices3d1d)(0,4) = 24; (*tIndices3d1d)(0,5) = 45;

// Subelement 2
(*tIndices3d1d)(0,6) = 3; (*tIndices3d1d)(0,7) = 13; (*tIndices3d1d)(0,8) = 16;
(*tIndices3d1d)(0,9) = 21; (*tIndices3d1d)(0,10) = 24; (*tIndices3d1d)(0,11) = 45;

// Subelement 3
(*tIndices3d1d)(0,12) = 12; (*tIndices3d1d)(0,13) = 13; (*tIndices3d1d)(0,14) = 14;
(*tIndices3d1d)(0,15) = 21; (*tIndices3d1d)(0,16) = 22; (*tIndices3d1d)(0,17) = 45;

// Subelement 4
(*tIndices3d1d)(0,18) = 7; (*tIndices3d1d)(0,19) = 14; (*tIndices3d1d)(0,20) = 15;
(*tIndices3d1d)(0,21) = 22; (*tIndices3d1d)(0,22) = 23; (*tIndices3d1d)(0,23) = 45;

// Subelement 5
(*tIndices3d1d)(0,24) = 11; (*tIndices3d1d)(0,25) = 15; (*tIndices3d1d)(0,26) = 16;
(*tIndices3d1d)(0,27) = 25; (*tIndices3d1d)(0,28) = 26; (*tIndices3d1d)(0,29) = 46;

// Subelement 6
(*tIndices3d1d)(0,30) = 2; (*tIndices3d1d)(0,31) = 16; (*tIndices3d1d)(0,32) = 18;
(*tIndices3d1d)(0,33) = 26; (*tIndices3d1d)(0,34) = 27; (*tIndices3d1d)(0,35) = 46;

// Subelement 7
(*tIndices3d1d)(0,36) = 10; (*tIndices3d1d)(0,37) = 17; (*tIndices3d1d)(0,38) = 18;
(*tIndices3d1d)(0,39) = 27; (*tIndices3d1d)(0,40) = 28; (*tIndices3d1d)(0,41) = 46;

// Subelement 8
(*tIndices3d1d)(0,42) = 6; (*tIndices3d1d)(0,43) = 15; (*tIndices3d1d)(0,44) = 17;
(*tIndices3d1d)(0,45) = 25; (*tIndices3d1d)(0,46) = 28; (*tIndices3d1d)(0,47) = 46;

// Subelement 9
(*tIndices3d1d)(0,48) = 10; (*tIndices3d1d)(0,49) = 17; (*tIndices3d1d)(0,50) = 18;
(*tIndices3d1d)(0,51) = 30; (*tIndices3d1d)(0,52) = 31; (*tIndices3d1d)(0,53) = 47;

// Subelement 10
(*tIndices3d1d)(0,54) = 1; (*tIndices3d1d)(0,55) = 18; (*tIndices3d1d)(0,56) = 19;
(*tIndices3d1d)(0,57) = 31; (*tIndices3d1d)(0,58) = 32; (*tIndices3d1d)(0,59) = 47;

// Subelement 11
(*tIndices3d1d)(0,60) = 9; (*tIndices3d1d)(0,61) = 19; (*tIndices3d1d)(0,62) = 20;
(*tIndices3d1d)(0,63) = 29; (*tIndices3d1d)(0,64) = 32; (*tIndices3d1d)(0,65) = 47;

// Subelement 12
(*tIndices3d1d)(0,66) = 5; (*tIndices3d1d)(0,67) = 17; (*tIndices3d1d)(0,68) = 20;
(*tIndices3d1d)(0,69) = 29; (*tIndices3d1d)(0,70) = 30; (*tIndices3d1d)(0,71) = 47;

// Subelement 13
(*tIndices3d1d)(0,72) = 12; (*tIndices3d1d)(0,73) = 13; (*tIndices3d1d)(0,74) = 14;
(*tIndices3d1d)(0,75) = 33; (*tIndices3d1d)(0,76) = 34; (*tIndices3d1d)(0,77) = 48;

// Subelement 14
(*tIndices3d1d)(0,78) = 4; (*tIndices3d1d)(0,79) = 13; (*tIndices3d1d)(0,80) = 19;
(*tIndices3d1d)(0,81) = 34; (*tIndices3d1d)(0,82) = 35; (*tIndices3d1d)(0,83) = 48;

// Subelement 15
(*tIndices3d1d)(0,84) = 9; (*tIndices3d1d)(0,85) = 19; (*tIndices3d1d)(0,86) = 20;
(*tIndices3d1d)(0,87) = 35; (*tIndices3d1d)(0,88) = 36; (*tIndices3d1d)(0,89) = 48;

// Subelement 16
(*tIndices3d1d)(0,90) = 8; (*tIndices3d1d)(0,91) = 14; (*tIndices3d1d)(0,92) = 20;
(*tIndices3d1d)(0,93) = 33; (*tIndices3d1d)(0,94) = 36; (*tIndices3d1d)(0,95) = 48;

// Subelement 17
(*tIndices3d1d)(0,96) = 7; (*tIndices3d1d)(0,97) = 14; (*tIndices3d1d)(0,98) = 15;
(*tIndices3d1d)(0,99) = 37; (*tIndices3d1d)(0,100) = 38; (*tIndices3d1d)(0,101) = 49;

// Subelement 18
(*tIndices3d1d)(0,102) = 6; (*tIndices3d1d)(0,103) = 15; (*tIndices3d1d)(0,104) = 17;
(*tIndices3d1d)(0,105) = 38; (*tIndices3d1d)(0,106) = 39; (*tIndices3d1d)(0,107) = 49;

// Subelement 19
(*tIndices3d1d)(0,108) = 5; (*tIndices3d1d)(0,109) = 17; (*tIndices3d1d)(0,110) = 20;
(*tIndices3d1d)(0,111) = 39; (*tIndices3d1d)(0,112) = 40; (*tIndices3d1d)(0,113) = 49;

// Subelement 20
(*tIndices3d1d)(0,114) = 8; (*tIndices3d1d)(0,115) = 14; (*tIndices3d1d)(0,116) = 20;
(*tIndices3d1d)(0,117) = 37; (*tIndices3d1d)(0,118) = 40; (*tIndices3d1d)(0,119) = 49;

// Subelement 21
(*tIndices3d1d)(0,120) = 3; (*tIndices3d1d)(0,121) = 13; (*tIndices3d1d)(0,122) = 16;
(*tIndices3d1d)(0,123) = 41; (*tIndices3d1d)(0,124) = 42; (*tIndices3d1d)(0,125) = 50;

// Subelement 22
(*tIndices3d1d)(0,126) = 2; (*tIndices3d1d)(0,127) = 16; (*tIndices3d1d)(0,128) = 18;
(*tIndices3d1d)(0,129) = 42; (*tIndices3d1d)(0,130) = 43; (*tIndices3d1d)(0,131) = 50;

// Subelement 23
(*tIndices3d1d)(0,132) = 1; (*tIndices3d1d)(0,133) = 18; (*tIndices3d1d)(0,134) = 19;
(*tIndices3d1d)(0,135) = 43; (*tIndices3d1d)(0,136) = 44; (*tIndices3d1d)(0,137) = 50;

// Subelement 24
(*tIndices3d1d)(0,138) = 4; (*tIndices3d1d)(0,139) = 13; (*tIndices3d1d)(0,140) = 19;
(*tIndices3d1d)(0,141) = 41; (*tIndices3d1d)(0,142) = 44; (*tIndices3d1d)(0,143) = 50;

//Offsets
(*tOffsets3d1d)(0,0) = 0; (*tOffsets3d1d)(0,1) = 6; (*tOffsets3d1d)(0,2) = 12; (*tOffsets3d1d)(0,3) = 18;
(*tOffsets3d1d)(0,4) = 24; (*tOffsets3d1d)(0,5) = 30; (*tOffsets3d1d)(0,6) = 36; (*tOffsets3d1d)(0,7) = 42;
(*tOffsets3d1d)(0,8) = 48; (*tOffsets3d1d)(0,9) = 54; (*tOffsets3d1d)(0,10) = 60; (*tOffsets3d1d)(0,11) = 66;
(*tOffsets3d1d)(0,12) = 72; (*tOffsets3d1d)(0,13) = 78; (*tOffsets3d1d)(0,14) = 84; (*tOffsets3d1d)(0,15) = 90;
(*tOffsets3d1d)(0,16) = 96; (*tOffsets3d1d)(0,17) = 102; (*tOffsets3d1d)(0,18) = 108; (*tOffsets3d1d)(0,19) = 114;
(*tOffsets3d1d)(0,20) = 120; (*tOffsets3d1d)(0,21) = 126; (*tOffsets3d1d)(0,22) = 132; (*tOffsets3d1d)(0,23) = 138;
(*tOffsets3d1d)(0,24) = 144;
