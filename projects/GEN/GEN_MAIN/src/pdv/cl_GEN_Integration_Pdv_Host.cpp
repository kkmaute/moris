#include "cl_GEN_Integration_Pdv_Host.hpp"
#include "cl_GEN_Pdv_Value.hpp"
#include "cl_GEN_Pdv_Intersection.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Integration_Pdv_Host::Integration_Pdv_Host(uint aNodeIndex,
                                                   const Matrix<DDRMat>& aCoordinates,
                                                   const Cell<PDV_Type>& aPdvTypes,
                                                   uint aStartingGlobalIndex)
        : mNodeIndex(aNodeIndex),
          mCoordinates(aCoordinates),
          mStartingGlobalIndex(aStartingGlobalIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Pdv_Host::~Integration_Pdv_Host()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Integration_Pdv_Host::create_pdv(PDV_Type aPdvType, real aPdvVal)
        {
            // Check PDV type
            MORIS_ASSERT(static_cast<uint>(aPdvType) < 3, 
                    "Tried to call Pdv_Host.create_pdv() using a fixed value with a PDV type that doesn't exist on this host");

            // Create a pdv with pdv value
            mPdvs(static_cast<uint>(aPdvType)) = std::make_shared< Pdv_Value >(aPdvVal);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Integration_Pdv_Host::create_pdv(PDV_Type aPdvType, GEN_Geometry_Object* aIntersection, uint aDimension)
        {
            // Create a pdv with an intersection pointer
            mPdvs(static_cast<uint>(aPdvType)) = std::make_shared<Pdv_Intersection>(aIntersection, aDimension);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Integration_Pdv_Host::is_active_type(PDV_Type aPdvType)
        {
            // return if active PDV
            return mPdvs(static_cast<uint>(aPdvType))->mIsActive;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Integration_Pdv_Host::get_global_index_for_pdv_type(PDV_Type aPdvType)
        {
            // Return id from map
            return mStartingGlobalIndex + static_cast<uint>(aPdvType);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Integration_Pdv_Host::get_starting_global_index()
        {
            return mStartingGlobalIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Integration_Pdv_Host::get_pdv_value(PDV_Type aPdvType)
        {
            // Return value
            return mPdvs(static_cast<uint>(aPdvType))->get_value(mNodeIndex, mCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Integration_Pdv_Host::get_all_sensitivities(Matrix<DDRMat>& aSensitivities)
        {
            aSensitivities.resize(0, aSensitivities.n_cols());
            if (mPdvs.size() > 0)
            {
                uint tActivePdvIndex = 0;
                Matrix<DDRMat> tSensitivities(0, 0);
                for (uint tPdvIndex = 0; tPdvIndex < mPdvs.size(); tPdvIndex++)
                {
                    if (mPdvs(tPdvIndex)->mIsActive)
                    {
                        mPdvs(tPdvIndex)->get_sensitivity(mNodeIndex, mCoordinates, tSensitivities);
                        aSensitivities.resize(tActivePdvIndex + 1, tSensitivities.length());
                        aSensitivities.set_row(tActivePdvIndex++, tSensitivities);
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
