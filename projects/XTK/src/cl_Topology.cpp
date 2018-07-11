#include "cl_Topology.hpp"

// ----------------------------------------------------------------------------
// Default constructor
moris::xtk::Topology::Topology()
{
}

moris::xtk::Topology::~Topology()
{
}

// ----------------------------------------------------------------------------
void
moris::xtk::Topology::set_topology(moris::lint                          aDim1,
                                   moris::lint                          aDim2,
                                   moris::xtk::Connectivity & aTopo)
{
    mTopologies( aDim1 )( aDim2 ) = aTopo;
}


// ----------------------------------------------------------------------------

moris::xtk::Connectivity*
moris::xtk::Topology::operator()(
        moris::lint const  d,
        moris::lint const  d_prime)
{
    return &mTopologies( d )( d_prime);
}


//moris::xtk::Connectivity*
//moris::xtk::Topology::operator()(
//        moris::lint const  d,
//        moris::lint const  d_prime) const
//{
//    return &mTopologies( d )( d_prime);
//}
// ----------------------------------------------------------------------------

//moris::xtk::Connectivity
//moris::xtk::Topology::operator()(
//        moris::lint const  d,
//        moris::lint const  d_prime)
//{
//
//    return mTopologies( d )( d_prime );
//
//}

// ------------------------------------------------------------------------------

