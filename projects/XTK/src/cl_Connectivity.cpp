#include "cl_Connectivity.hpp"

// ----------------------------------------------------------------------------
// Default constructor
moris::xtk::Connectivity::Connectivity()
{
}

// ----------------------------------------------------------------------------
// Default constructor
moris::xtk::Connectivity::~Connectivity()
{
}

// ----------------------------------------------------------------------------
//void
//moris::xtk::Connectivity::set_connectivity(moris::Cell< moris::Cell< moris::lint > > & aConn)
//{
//    mConn = aConn;
//}



// ----------------------------------------------------------------------------

void
moris::xtk::Connectivity::set_connectivity(moris::Mat<moris::uint> & aConn)
{
    mConn2 = aConn;
}

// ----------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::num_entities()
{
    return mConn.size();
}

// ----------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::num_entities() const
{
    return mConn.size();
}

// ----------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::num_conn( moris::lint & arow)
{
    return mConn(arow).size();
}

// ----------------------------------------------------------------------------

moris::uint
moris::xtk::Connectivity::num_conn( moris::uint & arow)
{
    return mConn(arow).size();
}

// ----------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::num_conn( moris::lint & arow) const
{
    return mConn(arow).size();
}

// ----------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::operator()(
        moris::lint const & arow,
        moris::lint const & acol)
{
    return mConn(arow)(acol);
}

// ----------------------------------------------------------------------------

moris::lint const
moris::xtk::Connectivity::operator()(
        moris::lint const & arow,
        moris::lint const & acol) const
{
    return mConn(arow)(acol);
}

// ------------------------------------------------------------------------------

moris::lint
moris::xtk::Connectivity::operator()(
        moris::uint const & arow,
        moris::uint const & acol)
{
    return mConn(arow)(acol);
}


moris::lint
moris::xtk::Connectivity::get_node(moris::uint    aRow,
                                   moris::uint    aCol)
{
    return mConn(aRow)(aCol);
}

moris::Cell<moris::uint>
moris::xtk::Connectivity::get_row(moris::uint  & aRow)
{
    return mConn(aRow);
}

//void
//test_func()
//{
//    moris::Cell<moris::Cell<moris::lint>> parent0Dto1D;
//
//        parent0Dto1D[0][0] =  1;  parent0Dto1D[0][1] =  4;
//        parent0Dto1D[1][0] =  1;  parent0Dto1D[1][1] =  2;  parent0Dto1D[1][2] =  5;
//        parent0Dto1D[2][0] =  2;  parent0Dto1D[2][1] =  3;  parent0Dto1D[2][2] =  7;  parent0Dto1D[2][3] = 10;
//        parent0Dto1D[3][0] =  4;  parent0Dto1D[3][1] = 12;
//        parent0Dto1D[4][0] =  5;  parent0Dto1D[4][1] =  6;
//        parent0Dto1D[5][0] =  6;  parent0Dto1D[5][1] =  7;  parent0Dto1D[5][2] =  8;
//        parent0Dto1D[6][0] =  8;  parent0Dto1D[6][1] =  9;
//        parent0Dto1D[7][0] =  9;  parent0Dto1D[7][1] = 10;  parent0Dto1D[7][2] = 11;
//        parent0Dto1D[8][0] = 11;  parent0Dto1D[8][1] = 12;
//
//        moris::xtk::Connectivity b0Dto1D;
//
//        b0Dto1D.set_connectivity( parent0Dto1D );
//}
