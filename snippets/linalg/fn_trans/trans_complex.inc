moris::Mat< moris::cplx > A( 2, 2 );

A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0};
A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0};

moris::Mat< moris::cplx > B( 2, 2 );
B = moris::trans( A );