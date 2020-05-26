#ifndef MORIS_CL_GEN_PDV_VALUE_HPP
#define MORIS_CL_GEN_PDV_VALUE_HPP

#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Value : public Pdv
        {

        private:
            real mValue;

        public:
            /**
             * Constructor
             *
             * @param aValue Constant value for this PDV
             */
            Pdv_Value(real aValue);

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            real get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * Get the PDV sensitivity with respect to ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @param aSensitivities Matrix of sensitivities to be returned
             */
            void get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_PDV_VALUE_HPP
