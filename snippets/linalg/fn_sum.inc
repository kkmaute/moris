moris::Mat< moris::real > aMatr (3, 3 );

aMatr( 0, 0 ) = 1.0; aMatr( 0, 1 ) = 2.0; aMatr( 0, 2 ) = 3.0;
aMatr( 1, 0 ) = 4.0; aMatr( 1, 1 ) = 5.0; aMatr( 1, 2 ) = 6.0;
aMatr( 2, 0 ) = 7.0; aMatr( 2, 1 ) = 8.0; aMatr( 2, 2 ) = 9.0;

moris::real sumr = moris::sum( aMatr );
