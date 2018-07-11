std::vector< std::vector< moris::lint > > parent0Dto1D( 9, std::vector< moris::lint >( 2 ));
//
parent0Dto1D[0][0] =  0;  parent0Dto1D[0][1] =  3;
parent0Dto1D[1][0] =  0;  parent0Dto1D[1][1] =  1;  parent0Dto1D[1].push_back( 4 );
parent0Dto1D[2][0] =  1;  parent0Dto1D[2][1] =  2;  parent0Dto1D[2].push_back( 6 );  parent0Dto1D[2].push_back( 9 );
parent0Dto1D[3][0] =  3;  parent0Dto1D[3][1] = 11;
parent0Dto1D[4][0] =  4;  parent0Dto1D[4][1] =  5;
parent0Dto1D[5][0] =  5;  parent0Dto1D[5][1] =  6;  parent0Dto1D[5].push_back( 7 );
parent0Dto1D[6][0] =  7;  parent0Dto1D[6][1] =  8;
parent0Dto1D[7][0] =  8;  parent0Dto1D[7][1] =  9;  parent0Dto1D[7].push_back( 10 );
parent0Dto1D[8][0] = 10;  parent0Dto1D[8][1] = 11;

moris::xtk::Connectivity b0Dto1D;

b0Dto1D.set_connectivity( parent0Dto1D );

moris::lint arow = 2;
moris::lint acol = 2;

moris::lint aValue1 = b0Dto1D.num_entities();
moris::lint aValue2 = b0Dto1D.num_conn( arow );
moris::lint aValue3 = b0Dto1D( arow, acol );
