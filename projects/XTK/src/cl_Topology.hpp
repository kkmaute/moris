#ifndef MORIS_XTK_TOPOLOGY_HPP_
#define MORIS_XTK_TOPOLOGY_HPP_

// MORIS library header files.


// MORIS project header files.
#include "core.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "assert.hpp"
#include "cl_Connectivity.hpp"

namespace moris
{
namespace xtk
{
    class Topology
    {
    public:

        /**
         * Constructor. This Class serves as a rough template for nodal cluster information for the enrichment strategy useage.
         * Some topological dimensions are not available!
         *
         * Examples of use:
         * @include /xtk/cl_xtk_Topology1.inc
         */
        Topology();

        /**
         * Destructor.
         */
        ~Topology();

        /**
         * set_topology operator for connectivity aDim1 -> aDim2.
         *
         * @param aDim1       Topological dimension aDim1.
         *
         * @param aDim2       Topological dimension aDim2.
         *
         * @param aTopo       Connectivity Topo for specific connectivity.
         *
         * Example:
         * @include XTK/src/fn_GetEnrichmentTable.inc
         *
         */
        void
        set_topology(moris::lint                aDim1,
                     moris::lint                aDim2,
                     moris::xtk::Connectivity & aTopo);




        /**
         * Constant accessor operator for connectivity d -> d_prime.
         *
         * @param d       Topological dimension d.
         * @param d_prime Topological dimension d_prime.
         * @return Mesh connectivity.
         */
        moris::xtk::Connectivity*
        operator()(
                moris::lint const  d,
                moris::lint const  d_prime);

        /**
         * Constant accessor operator for connectivity d -> d_prime.
         *
         * @param d       Topological dimension d.
         * @param d_prime Topological dimension d_prime.
         * @return Mesh connectivity.
         */
//        moris::xtk::Connectivity*
//        operator()(
//                moris::lint const  d,
//                moris::lint const  d_prime) const;

        /**
         * Constant accessor operator for connectivity d -> d_prime.
         *
         * @param d       Topological dimension d.
         * @param d_prime Topological dimension d_prime.
         * @return Mesh connectivity.
         */
//        moris::xtk::Connectivity
//        operator()(
//                moris::lint const  d,
//                moris::lint const  d_prime);


    private:

        moris::Cell< moris::Mat<moris::xtk::Connectivity > > mTopologies;

    };

}   // namespace mesh
}   // namespace moris

#endif  /* MORIS_XTK_TOPOLOGY_HPP_ */
