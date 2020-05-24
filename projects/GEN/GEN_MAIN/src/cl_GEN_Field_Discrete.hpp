//
// Created by christopherson on 5/19/20.
//

#ifndef MORIS_CL_GEN_FIELD_DISCRETE_HPP
#define MORIS_CL_GEN_FIELD_DISCRETE_HPP

#include "cl_GEN_Field_Base.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Discrete : virtual public Field
        {
        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Discrete();

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real evaluate_field_value(      uint            aIndex,
                                      const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aCoordinates vector of coordinate values
             * @return Field value
             */
            virtual real evaluate_field_value(uint aIndex) = 0;

        private:

            /**
             * Given a node coordinate, returns a matrix of all sensitivities
             *
             * @param aIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            void evaluate_all_sensitivities(      uint            aIndex,
                                            const Matrix<DDRMat>& aCoordinates,
                                                  Matrix<DDRMat>& aSensitivities);

            /**
             * Given a node index, returns a matrix of all sensitivities
             *
             * @param aCoordinates Vector of coordinate values
             * @param aSensitivity Matrix of sensitivities
             */
            virtual void evaluate_all_sensitivities(uint            aIndex,
                                                    Matrix<DDRMat>& aSensitivities) = 0;

        };
    }
}


#endif //MORIS_CL_GEN_FIELD_DISCRETE_HPP
