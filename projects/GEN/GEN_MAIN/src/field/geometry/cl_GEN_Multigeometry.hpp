#ifndef MORIS_CL_GEN_MULTIGEOMETRY_HPP
#define MORIS_CL_GEN_MULTIGEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {
        class Multigeometry : public Geometry
        {

        private:
            Cell<std::shared_ptr<Geometry>> mGeometries;

        public:

            /**
             * Multigeometry constructor
             *
             * @param aGeometries Created geometries
             */
            Multigeometry(Cell<std::shared_ptr<Geometry>> aGeometries);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real get_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index or coordinate, returns a matrix of all sensitivities
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Adds a geometry to this multigeometry.
             *
             * @param aGeometry Geometry to add
             */
            void add_geometry(std::shared_ptr<Geometry> aGeometry);

        };
    }
}

#endif //MORIS_CL_GEN_MULTIGEOMETRY_HPP
