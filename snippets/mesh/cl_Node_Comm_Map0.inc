// Create Node LIDs for elemental communication map.
moris::Mat<moris::size_t> nodes1 = {{0,5,1,4,2,3}};

// Create node communication map.
moris::mesh::Node_Comm_Map ncmap1(nodes1);
