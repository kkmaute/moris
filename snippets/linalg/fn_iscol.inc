moris::Mat< moris::real > A( 3, 3 );
moris::Mat< moris::real > B( 3, 1 );
moris::Mat< moris::real > C( 1, 3 );
moris::Mat< moris::real > D( 1, 1 );

A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

B( 0, 0 ) = 1.0; B( 1, 0 ) = 2.0; B( 2, 0 ) = 3.0;

C( 0, 0 ) = 1.0; C( 0, 1 ) = 2.0; C( 0, 2 ) = 3.0;

D( 0, 0 ) = 1.0;

bool e = moris::iscol( A );
bool f = moris::iscol( B );
bool g = moris::iscol( C );
bool h = moris::iscol( D );