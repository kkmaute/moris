moris::Mat< moris::real > w( 3,1 );
w( 0,0 ) = 1.0; w( 1,0 ) = 2.0; w( 2,0 ) = 3.0;

moris::Mat< moris::real > C = moris::diag( w );