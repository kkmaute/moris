

    moris::xtk::Topology assoctopo;

    assoctopo.set_topology(D1,D1,a1Dto1D);
    assoctopo.set_topology(D3,D3,a3Dto3D);

    moris::xtk::Topology childtopo;

    childtopo.set_topology(D1,D3,c1Dto3D);

    moris::Mat< moris::lint > EnrichmentTable;

    moris::Mat< moris::lint > ConnPair;

    moris::xtk::GetEnrichmentTable(EnrichmentTable, parenttopo, assoctopo, childtopo, IntLevels);
