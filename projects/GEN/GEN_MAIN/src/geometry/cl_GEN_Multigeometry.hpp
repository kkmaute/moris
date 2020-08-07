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
            std::string mID;

        public:

            /**
             * Multigeometry constructor
             *
             * @param aGeometries Already-created geometries
             * @param aID Name ID for this multigeometry
             */
            Multigeometry(Cell<std::shared_ptr<Geometry>> aGeometries,
                          std::string aID);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node index or coordinate, returns a matrix of relevant sensitivities
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            void evaluate_sensitivity(uint                  aNodeIndex,
                                      const Matrix<DDRMat>& aCoordinates,
                                      Matrix<DDRMat>&       aSensitivities);

            /**
             * Adds a geometry to this multigeometry.
             *
             * @param aGeometry Geometry to add
             */
            void add_geometry(std::shared_ptr<Geometry> aGeometry);

            /**
             * Gets the ID of this multigeometry.
             *
             * @return ID
             */
            std::string get_id();

        private:

            /**
             * Given a node coordinate @param aCoordinates, the function returns a matrix of sensitivities of the
             * geometry location with respect to the ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Matrix of sensitivities
             */
            void evaluate_all_sensitivities(uint                  aNodeIndex,
                                            const Matrix<DDRMat>& aCoordinates,
                                            Matrix<DDRMat>&       aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_MULTIGEOMETRY_HPP
