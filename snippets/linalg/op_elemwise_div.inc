moris::Mat< moris::real > Am( 2, 2 );
moris::Mat< moris::real > Bm( 2, 2 );

Am( 0, 0 ) = 9.0;   Am( 0, 1 ) = 12.0;
Am( 1, 0 ) = 400.0; Am( 1, 1 ) = 50.0;

Bm( 0, 0 ) = 3.0;   Bm( 0, 1 ) = 6.0;
Bm( 1, 0 ) = 100.0; Bm( 1, 1 ) = 10.0;

moris::Mat< moris::real > Cm( 2, 2 );
Cm = Am / Bm;

moris::Mat< moris::real > Asc( 2, 2 );

Asc( 0, 0 ) = 9.0;   Asc( 0, 1 ) = 12.0;
Asc( 1, 0 ) = 300.0; Asc( 1, 1 ) = 30.0;

moris::real Bsc = 3.0;

moris::Mat< moris::real > Csc( 2, 2 );
Csc = Asc / Bsc;