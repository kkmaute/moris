// GEN
#include "cl_GEN_Geometric_Proximity.hpp"

namespace moris
{
    namespace ge
    {
        Geometric_Proximity::Geometric_Proximity() {}

        Geometric_Proximity::Geometric_Proximity( moris_index const & aNumGeometries )
                : mAssociatedVertexIndex( MORIS_INDEX_MAX )
                , mGeometricProximity( 1, aNumGeometries, MORIS_INDEX_MAX )
        {
        }

        //-----------------------------------------------------------------------------------------

        Geometric_Proximity::~Geometric_Proximity() {}

        //-----------------------------------------------------------------------------------------

        void
        Geometric_Proximity::set_geometric_proximity( moris_index aGeometricProximity, moris_index aGeometryIndex )
        {
            MORIS_ASSERT( aGeometryIndex < (moris_index)mGeometricProximity.numel(), "Geometry index out of bounds" );
            mGeometricProximity( aGeometryIndex ) = aGeometricProximity;
        }

        //-----------------------------------------------------------------------------------------

        moris_index
        Geometric_Proximity::get_geometric_proximity( moris_index aGeometryIndex )
        {
            return mGeometricProximity( aGeometryIndex );
        }

        //-----------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
