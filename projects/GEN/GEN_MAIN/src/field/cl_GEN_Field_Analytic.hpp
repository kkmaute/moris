#ifndef MORIS_CL_GEN_FIELD_ANALYTIC_HPP
#define MORIS_CL_GEN_FIELD_ANALYTIC_HPP

#include "cl_GEN_Field.hpp"

namespace moris
{
    namespace ge
    {
        class Field_Analytic : virtual public Field
        {
        public:

            /**
             * Trivial constructor, necessary for clean virtual inheritance without default constructor in base class
             */
            Field_Analytic();

            /**
             * Given a node index or coordinate, returns the field value.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Field value
             */
            real get_field_value(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns the field value
             *
             * @param aCoordinates vector of coordinate values
             * @return Field value
             */
            virtual real get_field_value(const Matrix<DDRMat>& aCoordinates) = 0;

            /**
             * Given a node index or coordinate, returns a matrix of all sensitivities.
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            const Matrix<DDRMat>& get_field_sensitivities(
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Given a node coordinate, returns a matrix of all sensitivities
             *
             * @param aCoordinates Vector of coordinate values
             * @return Matrix of sensitivities
             */
            virtual const Matrix<DDRMat>& get_field_sensitivities(const Matrix<DDRMat>& aCoordinates) = 0;

        };
    }
}

#endif //MORIS_CL_GEN_FIELD_ANALYTIC_HPP
