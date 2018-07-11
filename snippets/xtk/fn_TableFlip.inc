moris::Mat< moris::lint > ConTable( 17, 4 );

ConTable( 0, 0 ) =  38; ConTable( 0, 1 ) =  37; ConTable( 0, 2 ) =  -1; ConTable( 0, 3 ) = -1;
ConTable( 1, 0 ) =  13; ConTable( 1, 1 ) =  15; ConTable( 1, 2 ) =  -1; ConTable( 1, 3 ) = -1;
ConTable( 2, 0 ) =   0; ConTable( 2, 1 ) =   1; ConTable( 2, 2 ) =  -1; ConTable( 2, 3 ) = -1;
ConTable( 3, 0 ) =   0; ConTable( 3, 1 ) =  -1; ConTable( 3, 2 ) =  -1; ConTable( 3, 3 ) =  1;
ConTable( 4, 0 ) =  26; ConTable( 4, 1 ) =  -1; ConTable( 4, 2 ) =  -1; ConTable( 4, 3 ) = 25;
ConTable( 5, 0 ) =  37; ConTable( 5, 1 ) =  -1; ConTable( 5, 2 ) =  -1; ConTable( 5, 3 ) = 37;
ConTable( 6, 0 ) =  -1; ConTable( 6, 1 ) =  38; ConTable( 6, 2 ) =  37; ConTable( 6, 3 ) = -1;
ConTable( 7, 0 ) =  -1; ConTable( 7, 1 ) =  25; ConTable( 7, 2 ) =  25; ConTable( 7, 3 ) = -1;
ConTable( 8, 0 ) =  -1; ConTable( 8, 1 ) =   1; ConTable( 8, 2 ) =   0; ConTable( 8, 3) =  -1;
ConTable( 9, 0 ) =  -1; ConTable( 9, 1 ) =  -1; ConTable( 9, 2 ) =   0; ConTable( 9, 3) =   1;
ConTable( 10, 0 ) =  1; ConTable( 10, 1 ) = -1; ConTable( 10, 2 ) = -1; ConTable( 10, 3) = -1;
ConTable( 11, 0 ) = 25; ConTable( 11, 1 ) = -1; ConTable( 11, 2 ) = -1; ConTable( 11, 3) = -1;
ConTable( 12, 0 ) = -1; ConTable( 12, 1 ) = 14; ConTable( 12, 2 ) = -1; ConTable( 12, 3) = -1;
ConTable( 13, 0 ) = -1; ConTable( 13, 1 ) =  0; ConTable( 13, 2 ) = -1; ConTable( 13, 3) = -1;
ConTable( 14, 0 ) = -1; ConTable( 14, 1 ) = 13; ConTable( 14, 2 ) = -1; ConTable( 14, 3) = -1;
ConTable( 15, 0 ) = -1; ConTable( 15, 1 ) = -1; ConTable( 15, 2 ) = 13; ConTable( 15, 3) = -1;
ConTable( 16, 0 ) = -1; ConTable( 16, 1 ) = -1; ConTable( 16, 2 ) = -1; ConTable( 16, 3) = 13;

moris::Mat< moris::lint > ConPair( 17, 1 );

ConPair( 0 ) =  1; ConPair( 1 ) =  1; ConPair( 2 ) =  1; ConPair( 3 ) =  1; ConPair( 4 ) =  1;
ConPair( 5 ) =  1; ConPair( 6 ) =  1; ConPair( 7 ) =  1; ConPair( 8 ) =  1; ConPair( 9 ) =  1;
ConPair( 10 ) = 0; ConPair( 11 ) = 0; ConPair( 12 ) = 0; ConPair( 13 ) = 0; ConPair( 14 ) = 0;
ConPair( 15 ) = 0; ConPair( 16 ) = 0;

moris::Mat< moris::lint > NewTable( 17, 4 );

NewTable.fill(-1);

moris::xtk::TableFlip( NewTable, ConTable, ConPair);
