/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Plane.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Analytic.hpp"

namespace moris::gen
{
    class Plane : public Field_Analytic< 2 > // FIXME plane should by default be 3D
    {
    private:
        real (Plane::*m_eval_field)(const Matrix<DDRMat>&) = nullptr;
        const Matrix<DDRMat>& (Plane::*m_eval_sensitivity)(const Matrix<DDRMat>&) = nullptr;

    public:
        /**
         * Constructor, sets the pointers to advs and constant parameters for evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         */
        Plane( ADV_ARG_TYPES )
               : Field_Analytic< 2 >( ADV_ARGS )
        {
            if ( mADVManager.get_determining_adv_ids().length() == 4 )
            {
                m_eval_field = &Plane::eval_field_2d;
                m_eval_sensitivity = &Plane::eval_sensitivity_2d;
            }
            else if ( mADVManager.get_determining_adv_ids().length() == 6 )
            {
                m_eval_field = &Plane::eval_field_3d;
                m_eval_sensitivity = &Plane::eval_sensitivity_3d;
            }
            else
            {
                MORIS_ERROR(false, "Incorrect number of parameters passed for construction of a GEN Plane.");
            }
        }

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
        Plane( real aXCenter,
               real aYCenter,
               real aZCenter,
               real aXNormal,
               real aYNormal,
               real aZNormal );

        /**
         * Constructor with only constant parameters, 2D
         *
         * @param aXCenter x-coordinate of the center of the plane
         * @param aYCenter y-coordinate of the center of the plane
         * @param aXNormal x normal for the plane
         * @param aYNormal y normal for the plane
         */
        Plane( real aXCenter,
               real aYCenter,
               real aXNormal,
               real aYNormal );

        /**
         * Given a node coordinate, returns the field value.
         *
         * @param aCoordinates Coordinate values
         * @return Field value
         */
        real get_field_value(const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node coordinate, evaluates the sensitivity of the field with respect to all of the
         * field variables.
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

    private:

        /**
         * 2D evaluation for get_voxel_id
         */
        real eval_field_2d(const Matrix<DDRMat>& aCoordinates);

        /**
         * 3D evaluation for get_voxel_id
         */
        real eval_field_3d(const Matrix<DDRMat>& aCoordinates);

        /**
         * 2D evaluation for get_dfield_dadvs
         */
        const Matrix<DDRMat>& eval_sensitivity_2d(const Matrix<DDRMat>& aCoordinates);

        /**
         * 3D evaluation for get_dfield_dadvs
         */
        const Matrix<DDRMat>& eval_sensitivity_3d(const Matrix<DDRMat>& aCoordinates);
    };
}
