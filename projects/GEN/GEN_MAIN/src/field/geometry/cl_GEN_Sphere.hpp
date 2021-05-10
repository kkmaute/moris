#ifndef MORIS_CL_GEN_SPHERE_HPP_
#define MORIS_CL_GEN_SPHERE_HPP_

#include <cmath>

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Sphere : public Geometry, public Field_Analytic
        {
        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations.
             *
             * @tparam Vector_Type Type of vector where ADVs are stored
             * @param aADVs ADV vector
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            template <typename Vector_Type>
            Sphere(Vector_Type&              aADVs,
                   Matrix<DDUMat>            aGeometryVariableIndices,
                   Matrix<DDUMat>            aADVIndices,
                   Matrix<DDRMat>            aConstants,
                   Geometry_Field_Parameters aParameters = {})
                    : Field(aADVs, aGeometryVariableIndices, aADVIndices, aConstants, aParameters)
                    , Geometry(aParameters)
            {
                MORIS_ERROR(aGeometryVariableIndices.length() + aConstants.length() == 4,
                            "A GEN Sphere must be created with a total of exactly 4 variables (ADVs + constants)");
            }

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the sphere
             * @param aYCenter y-coordiante of the center of the sphere
             * @param aZCenter z-coordinate of the center of the sphere
             * @param aRadius radius of the sphere
             * @param aParameters Additional parameters
             */
            Sphere(real                      aXCenter,
                   real                      aYCenter,
                   real                      aZCenter,
                   real                      aRadius,
                   Geometry_Field_Parameters aParameters = {});

            /**
             * Given a node coordinate, returns the field value.
             *
             * @param aCoordinates Coordinate values
             * @return Distance to this geometry
             */
            real get_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables.
             *
             * @param aCoordinates Coordinate values
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_dfield_dadvs(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given nodal coordinates, returns a vector of the field derivatives with respect to the nodal
             * coordinates.
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
             */
            void get_dfield_dcoordinates(
                    const Matrix<DDRMat>& aCoordinates,
                    Matrix<DDRMat>&       aSensitivities);

        };
    }
}

#endif /* MORIS_CL_GEN_SPHERE_HPP_ */
