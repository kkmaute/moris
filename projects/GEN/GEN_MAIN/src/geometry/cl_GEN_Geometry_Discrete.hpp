//
// Created by christopherson on 4/16/20.
//

#ifndef MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
#define MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry_Discrete
        {
        protected:
            Cell<real*> mGeometryVariables;
            Matrix<DDRMat> mConstantParameters;

            /**
             * Constructor, sets the pointers to advs and constant parameters for evaluations
             *
             * @param aADVs Reference to the full advs
             * @param aGeometryVariableIndices Indices of geometry variables to be filled by the ADVs
             * @param aADVIndices The indices of the ADV vector to fill in the geometry variables
             * @param aConstantParameters The constant parameters not filled by ADVs
             */
            Geometry_Discrete(Matrix<DDRMat>& aADVs, Matrix<DDUMat> aGeometryVariableIndices, Matrix<DDUMat> aADVIndices, Matrix<DDRMat> aConstantParameters);

            /**
             * Constructor for only constant parameters
             *
             * @param aConstantParameters The parameters that define this geometry
             */
            Geometry_Discrete(Matrix<DDRMat> aConstantParameters);

        public:
            /**
             * Destructor
             */
            ~Geometry_Discrete()
            {
            }

            /**
             * Given an index, the discrete geometry needs to return a field value.
             *
             * @param aEntityIndex the index of the field value
             * @return field value at the specified index
             */
            virtual real evaluate_field_value(moris_index aEntityIndex) = 0;

            /**
             * Given an index, the discrete geometry needs to return sensitivites with respect to the field value
             *
             * @param aEntityIndex the index of the field value
             * @return matrix of sensitivities
             */
            virtual Matrix<DDRMat> evaluate_sensitivity(moris_index aEntityIndex) = 0;

        };
    }
}

#endif //MORIS_CL_GEN_GEOMETRY_DISCRETE_HPP
