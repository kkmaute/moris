moris::Mat< moris::real > Am( 3, 3 );
moris::Mat< moris::real > Bm( 3, 3 );

Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 2.0; Am( 0, 2 ) = 3.0;
Am( 1, 0 ) = 4.0; Am( 1, 1 ) = 5.0; Am( 1, 2 ) = 6.0;
Am( 2, 0 ) = 7.0; Am( 2, 1 ) = 7.0; Am( 2, 2 ) = 9.0;

Bm( 0, 0 ) = 3.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 1.0;
Bm( 1, 0 ) = 6.0; Bm( 1, 1 ) = 5.0; Bm( 1, 2 ) = 4.0;
Bm( 2, 0 ) = 9.0; Bm( 2, 1 ) = 8.0; Bm( 2, 2 ) = 7.0;

moris::Mat< moris::uint > Cm = ( Am < Bm );

moris::Mat< moris::real > Ac( 3,1 );
moris::Mat< moris::real > Bc( 3,1 );

Ac( 0,0 ) = 1.0;
Ac( 1,0 ) = 0.0;
Ac( 2,0 ) = 3.0;

Bc( 0,0 ) = 1.0;
Bc( 1,0 ) = 2.0;
Bc( 2,0 ) = 1.0;

moris::Mat< moris::uint > Cc = ( Ac < Bc );