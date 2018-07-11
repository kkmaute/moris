moris::Mat< moris::real > A( 3, 3 );
moris::Mat< moris::real > L;

A( 0, 0 ) =  2.0; A( 0, 1 ) = -1.0; A( 0, 2 ) =  0.0;
A( 1, 0 ) = -1.0; A( 1, 1 ) =  2.0; A( 1, 2 ) = -1.0;
A( 2, 0 ) =  0.0; A( 2, 1 ) = -1.0; A( 2, 2 ) =  2.0;

moris::choll( L, A );