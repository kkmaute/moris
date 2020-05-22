#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv::Pdv(std::shared_ptr< GEN_Field > aFieldPointer,
            moris_index                  aEntityIndex )
        {
            // assign pdv value from the property pointer
            mValue = aFieldPointer->get_field_val_at_vertex( aEntityIndex );

            // flag this PDV_Type so the interface knows it is from a Field and therefore not changing
            mIsActive = false;
        }

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
