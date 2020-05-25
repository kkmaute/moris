#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv::Pdv(moris::real aPdvVal )
        {
            // assign pdv value directly
            mValue = aPdvVal;
        }

        //--------------------------------------------------------------------------------------------------------------

        Pdv::~Pdv()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv::get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return mValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv::get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

    }   // end ge namespace
}   // end moris namespace
