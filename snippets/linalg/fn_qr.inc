moris::Mat< moris::real > A( 4, 3 );
moris::Mat< moris::real > Q;
moris::Mat< moris::real > R;

A( 0, 0 ) = 1.0; A( 0, 1 ) = -1.0; A( 0, 2 ) =  4.0;
A( 1, 0 ) = 1.0; A( 1, 1 ) =  4.0; A( 1, 2 ) = -2.0;
A( 2, 0 ) = 1.0; A( 2, 1 ) =  4.0; A( 2, 2 ) =  2.0;
A( 3, 0 ) = 1.0; A( 3, 1 ) = -1.0; A( 3, 2 ) =  0.0;

moris::qr( Q, R, A);
