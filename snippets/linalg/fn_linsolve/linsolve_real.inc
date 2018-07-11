moris::Mat< moris::real > A( 3, 3 );
moris::Mat< moris::real > B( 3, 1 );

A( 0, 0 ) = 3.0; A( 0, 1 ) = 0.0; A( 0, 2 ) = 3.0;
A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
A( 2, 0 ) = 9.0; A( 2, 1 ) = 1.0; A( 2, 2 ) = 9.0;

B( 0, 0 ) = 0.0;
B( 1, 0 ) = 1.0;
B( 2, 0 ) = 2.0;

moris::Mat< moris::real > C = moris::solve( A,B );