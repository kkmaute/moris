#ifndef MORIS_CL_GEN_SUPERELLIPSOID_HPP
#define MORIS_CL_GEN_SUPERELLIPSOID_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {
        class Superellipsoid : public Geometry, public Field_Analytic
        {
        private:
            real mEpsilon = 1E-8;

        public:

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstants The constant field variables not filled by ADVs
             * @param aParameters Additional parameters
             */
            Superellipsoid(
                    Matrix<DDRMat>&  aADVs,
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
            Superellipsoid(
                    sol::Dist_Vector* aOwnedADVs,
                    Matrix<DDUMat>    aGeometryVariableIndices,
                    Matrix<DDUMat>    aADVIndices,
                    Matrix<DDRMat>    aConstants,
                    Field_Parameters  aParameters = {});

            /**
             * Constructor with only constant parameters
             *
             * @param aXCenter x-coordinate of the center of the superellipsoid
             * @param aYCenter y-coordiante of the center of the superellipsoid
             * @param aZCenter z-coordinate of the center of the superellipsoid
             * @param aXSemidiameter Superellipsoid semi-diameter in the x direction
             * @param aYSemidiameter Superellipsoid semi-diameter in the y direction
             * @param aZSemidiameter Superellipsoid semi-diameter in the z direction
             * @param aExponent Superellipsoid exponent
             * @param aParameters Additional parameters
             */
            Superellipsoid(
                    real             aXCenter,
                    real             aYCenter,
                    real             aZCenter,
                    real             aXSemidiameter,
                    real             aYSemidiameter,
                    real             aZSemidiameter,
                    real             aExponent,
                    Field_Parameters aParameters = {});

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

#endif //MORIS_CL_GEN_SUPERELLIPSOID_HPP
