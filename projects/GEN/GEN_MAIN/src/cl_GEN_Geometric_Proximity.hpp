#ifndef MORIS_CL_Geometric_Proximity_HPP_
#define MORIS_CL_Geometric_Proximity_HPP_

namespace moris
{
    namespace ge
    {
        class Geometric_Proximity
        {
            public:
            Geometric_Proximity(){};
            Geometric_Proximity(moris_index const & aNumGeometries)
            :mAssociatedVertexIndex(MORIS_INDEX_MAX),
             mGeometricProximity(1,aNumGeometries,MORIS_INDEX_MAX)
            {};
            ~Geometric_Proximity(){};

            void
            set_geometric_proximity(moris_index aGeometricProximity, moris_index aGeometryIndex)
            {
                MORIS_ASSERT(aGeometryIndex < mGeometricProximity.numel(), "Geometry index out of bounds");
                mGeometricProximity(aGeometryIndex) = aGeometricProximity;
            }

            moris_index
            get_geometric_proximity(moris_index aGeometryIndex)
            {
                return mGeometricProximity(aGeometryIndex);
            }

            moris_index mAssociatedVertexIndex;

            // Keeps track of a vertex proximity to each geometry ( NumVerts x NumGeometries)
            // 0 - G(x) < threshold
            // 1 - G(x) == threshold
            // 2 - G(x) > threshold
            // Max not set
            Matrix<IndexMat> mGeometricProximity;

        };
    }
}

#endif //MORIS_CL_Geometric_Proximity_HPP_

