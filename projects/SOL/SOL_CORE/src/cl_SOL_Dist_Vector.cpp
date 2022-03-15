#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"

namespace moris
{
    namespace sol
    {
        //--------------------------------------------------------------------------------------------------------------

        Dist_Vector::Dist_Vector(
                Dist_Map* aMapClass,
                bool      aManageMap )
                : mMap( aMapClass )
                , mManageMap( aManageMap )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Dist_Vector::~Dist_Vector()
        {
            if ( mManageMap )
            {
                delete mMap;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        sol::Dist_Map*
        Dist_Vector::Dist_Vector::get_map()
        {
            return mMap;
        }

        //--------------------------------------------------------------------------------------------------------------
    }    // namespace sol
}    // namespace moris
