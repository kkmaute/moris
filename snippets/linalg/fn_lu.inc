moris::Mat< moris::real > A( 3, 3 );
moris::Mat< moris::real > L;
moris::Mat< moris::real > U;
moris::Mat< moris::real > P;

A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 4.0;
A( 1, 0 ) = 3.0; A( 1, 1 ) = 8.0; A( 1, 2 ) = 14.0;
A( 2, 0 ) = 2.0; A( 2, 1 ) = 6.0; A( 2, 2 ) = 13.0;

moris::lu( L, U, P, A );