moris::Mat< moris::real > aMatr (3, 3 );

aMatr( 0, 0 ) = 2.33; aMatr( 0, 1 ) = 0.00; aMatr( 0, 2 ) = 0.00;
aMatr( 1, 0 ) = 0.00; aMatr( 1, 1 ) = 2.33; aMatr( 1, 2 ) = 0.00;
aMatr( 2, 0 ) = 0.00; aMatr( 2, 1 ) = 0.00; aMatr( 2, 2 ) = 2.33;

moris::real detr = moris::det( aMatr );