moris::Mat< moris::real > A( 3, 3 );
moris::Mat< moris::real > B( 3, 3 );

A( 0, 0 ) = 1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) = 1.0;
A( 1, 0 ) = 1.0; A( 1, 1 ) = 1.0; A( 1, 2 ) = 1.0;
A( 2, 0 ) = 1.0; A( 2, 1 ) = 1.0; A( 2, 2 ) = 1.0;

B( 0, 0 ) = 1.0; B( 0, 1 ) = 2.0; B( 0, 2 ) = 2.0;
B( 1, 0 ) = 2.0; B( 1, 1 ) = 1.0; B( 1, 2 ) = 2.0;
B( 2, 0 ) = 2.0; B( 2, 1 ) = 2.0; B( 2, 2 ) = 1.0;

moris::Mat< moris::uint > C = ( A == B );