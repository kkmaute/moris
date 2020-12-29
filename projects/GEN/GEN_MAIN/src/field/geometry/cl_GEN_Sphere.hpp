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
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Sphere(Matrix<DDRMat>&  aADVs,
                   Matrix<DDUMat>   aGeometryVariableIndices,
                   Matrix<DDUMat>   aADVIndices,
                   Matrix<DDRMat>   aConstants,
                   Field_Parameters aParameters = {});

            /**
             * Constructor, sets the field variable pointers to ADVs and constant parameters for evaluations.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Sphere(sol::Dist_Vector* aOwnedADVs,
                   Matrix<DDUMat>    aGeometryVariableIndices,
                   Matrix<DDUMat>    aADVIndices,
                   Matrix<DDRMat>    aConstants,
                   Field_Parameters  aParameters = {});

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the sphere
             * @param aYCenter y-coordiante of the center of the sphere
             * @param aZCenter z-coordinate of the center of the sphere
             * @param aRadius radius of the sphere
             * @param aParameters Additional parameters
             */
            Sphere(real aXCenter, real aYCenter, real aZCenter, real aRadius, Field_Parameters aParameters = {});

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
            const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates);

        };
    }
}

#endif /* MORIS_CL_GEN_SPHERE_HPP_ */
