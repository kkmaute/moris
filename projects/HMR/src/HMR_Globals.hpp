#ifndef SRC_HMR_HMR_GLOBALS_HPP_
#define SRC_HMR_HMR_GLOBALS_HPP_
#include "typedefs.hpp" //COR/src
#include "cl_Param_List.hpp"       //CON/src

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        //datatype for hmr paramater list
        typedef Param_List< boost::variant< sint, real, std::string  > > ParameterList;

        // global parameter defining max supported level.
        // Needs to be known during compile time.
        const uint   gMaxNumberOfLevels   = 21;

        // Per default, padding elements are not owned by any proc.
        // This constant can be changed, eg to 0, if desired.
        const uint   gNoProcOwner      = MORIS_UINT_MAX;

        // value to be used if a proc has no neighbor at the
        // specified position
        // Must be identical to value used in create_proc_cart()
        const uint   gNoProcNeighbor     = MORIS_UINT_MAX;

        // value to be used if an element has no child index
        // const uint   gNoChildIndex        = MORIS_UINT_MAX;

        // value to be used if no basis id was given
        const luint  gNoEntityID          = MORIS_LUINT_MAX;

        const uint   gBitsetSize          =  3*gMaxNumberOfLevels;

        const real   gEpsilon             = 1e-6;

        const uint   gNumberOfPatterns    = 5;



// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
#endif
