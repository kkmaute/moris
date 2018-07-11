moris::size_t D3 = 3;
moris::size_t D2 = 2;
moris::size_t D1 = 1;
moris::size_t D0 = 0;

std::vector< std::vector< moris::lint > > parent1Dto3D( 12, std::vector< moris::lint >( 1 ));
//
parent1Dto3D[0][0] =  0;
parent1Dto3D[1][0] =  0;  parent1Dto3D[1].push_back( 1 );
parent1Dto3D[2][0] =  0;  parent1Dto3D[2].push_back( 3 );
parent1Dto3D[3][0] =  0;
parent1Dto3D[4][0] =  1;
parent1Dto3D[5][0] =  1;
parent1Dto3D[6][0] =  1;  parent1Dto3D[6].push_back( 2 );
parent1Dto3D[7][0] =  2;
parent1Dto3D[8][0] =  2;
parent1Dto3D[9][0] =  2;  parent1Dto3D[9].push_back( 3 );
parent1Dto3D[10][0] = 3;
parent1Dto3D[11][0] = 3;

moris::xtk::Connectivity b1Dto3D;

b1Dto3D.set_connectivity( parent1Dto3D );

moris::xtk::Topology parenttopo;

parenttopo.set_topology(D0,D1,b0Dto1D);
parenttopo.set_topology(D1,D3,b1Dto3D);


moris::xtk::Connectivity TestConn = parenttopo( D0, D1 );
