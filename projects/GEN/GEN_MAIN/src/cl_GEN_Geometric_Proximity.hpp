#ifndef MORIS_CL_Geometric_Proximity_HPP_
#define MORIS_CL_Geometric_Proximity_HPP_

#include "cl_Matrix.hpp"
#include "typedefs.hpp"
namespace moris
{
    namespace ge
    {
        class Geometric_Proximity
        {
          public:
            //-----------------------------------------------------------------------------------------

            Geometric_Proximity();

            //-----------------------------------------------------------------------------------------

            Geometric_Proximity( moris_index const & aNumGeometries );

            //-----------------------------------------------------------------------------------------

            ~Geometric_Proximity();

            //-----------------------------------------------------------------------------------------

            void
            set_geometric_proximity( moris_index aGeometricProximity, moris_index aGeometryIndex );

            //-----------------------------------------------------------------------------------------

            moris_index
            get_geometric_proximity( moris_index aGeometryIndex );

            //-----------------------------------------------------------------------------------------

            moris_index mAssociatedVertexIndex = MORIS_INDEX_MAX;

            //-----------------------------------------------------------------------------------------

            // Keeps track of a vertex proximity to each geometry ( NumVerts x NumGeometries)
            // 0 - G(x) < threshold
            // 1 - G(x) == threshold
            // 2 - G(x) > threshold
            // Max not set
            Matrix< IndexMat > mGeometricProximity;
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_Geometric_Proximity_HPP_
