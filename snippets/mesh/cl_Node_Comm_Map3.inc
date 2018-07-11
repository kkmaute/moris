// Empty nodal communication map.
moris::Mat<moris::uint> nodes2;

//trying to make a node communications map out of an empty matrix results in an error
moris::mesh::Node_Comm_Map ncmap2(nodes2);
