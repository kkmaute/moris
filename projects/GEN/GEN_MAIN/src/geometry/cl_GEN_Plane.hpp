#ifndef MORIS_CL_GEN_PLANE_HPP
#define MORIS_CL_GEN_PLANE_HPP

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Plane : public Geometry_Analytic
        {
        private:
            real ( Plane:: * m_eval_field )(Matrix<DDRMat> const &) = nullptr;
            Matrix<DDRMat> ( Plane:: * m_eval_sensitivity )(Matrix<DDRMat> const &) = nullptr;

        public:
            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Plane(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor with only constant parameters, 3D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aZCenter z-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             * @param aZNormal z normal for the plane
             */
            Plane(real aXCenter, real aYCenter, real aZCenter, real aXNormal, real aYNormal, real aZNormal);

            /**
             * Constructor with only constant parameters, 2D
             *
             * @param aXCenter x-coordinate of the center of the plane
             * @param aYCenter y-coordinate of the center of the plane
             * @param aXNormal x normal for the plane
             * @param aYNormal y normal for the plane
             */
            Plane(real aXCenter, real aYCenter, real aXNormal, real aYNormal);

            /**
             * Given a node coordinate, the geometry needs to return the distance to the nearest function.
             *
             * @param aCoordinates vector of coordinate values
             * @return distance to nearest function
             */
            real evaluate_field_value(const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
             * Where each row represents a design variable and each column is x, y, z sensitivities
             *
             * @param aCoordinates vector of coordinate values
             * @return matrix of sensitivities
             */
            Matrix<DDRMat> evaluate_sensitivity(const Matrix<DDRMat>& aCoordinates);


        private:
            /**
             * 2D evaluation for evaluate_field_value
             */
            moris::real eval_field_2d(Matrix<DDRMat> const & aCoordinates);

            /**
             * 3D evaluation for evaluate_field_value
             */
            moris::real eval_field_3d(Matrix<DDRMat> const & aCoordinates);

            /**
             * 2D evaluation for evaluate_sensitivity
             */
            Matrix<DDRMat> eval_sensitivity_2d(Matrix<DDRMat> const & aCoordinates);

            /**
             * 3D evaluation for evaluate_sensitivity
             */
            Matrix<DDRMat> eval_sensitivity_3d(Matrix<DDRMat> const & aCoordinates);
        };
    }
}


#endif /* MORIS_CL_GEN_PLANE_HPP */
