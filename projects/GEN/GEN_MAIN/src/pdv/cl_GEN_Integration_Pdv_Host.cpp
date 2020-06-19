#include "cl_GEN_Integration_Pdv_Host.hpp"
#include "cl_GEN_Geometry.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Integration_Pdv_Host::Integration_Pdv_Host(uint aNodeIndex,
                                                   const Matrix<DDRMat>& aCoordinates,
                                                   uint aStartingGlobalIndex,
                                                   std::shared_ptr<Intersection> aIntersectionInfo)
        : mNodeIndex(aNodeIndex),
          mCoordinates(aCoordinates),
          mStartingGlobalIndex(aStartingGlobalIndex),
          mIntersection(aIntersectionInfo)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Pdv_Host::~Integration_Pdv_Host()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Integration_Pdv_Host::set_as_intersection(std::shared_ptr<Intersection> aIntersection)
        {
            mIntersection = aIntersection;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Integration_Pdv_Host::is_active_type(PDV_Type aPdvType)
        {
            // If there is intersection info, PDVs are active
            return (mIntersection != nullptr);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Integration_Pdv_Host::get_global_index_for_pdv_type(PDV_Type aPdvType)
        {
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
            return mCoordinates(static_cast<uint>(aPdvType));
        }

        //--------------------------------------------------------------------------------------------------------------

        void Integration_Pdv_Host::get_all_sensitivities(Matrix<DDRMat>& aSensitivities)
        {
            if (mIntersection != nullptr)
            {
                MORIS_ERROR(mIntersection->mGeometry->sensitivities_available(),
                            "Sensitivities of an intersection not implemented if a geometry does not have its sensitivities available.");

                // Geometry field sensitivities
                Matrix<DDRMat> tFieldSensitivityNode1(0, 0);
                Matrix<DDRMat> tFieldSensitivityNode2(0, 0);
                mIntersection->mGeometry->evaluate_sensitivity(mIntersection->mNodes(0)->mNodeIndex,
                                                                   mIntersection->mNodes(0)->mCoordinates,
                                                                   tFieldSensitivityNode1);
                mIntersection->mGeometry->evaluate_sensitivity(mIntersection->mNodes(1)->mNodeIndex,
                                                                   mIntersection->mNodes(1)->mCoordinates,
                                                                   tFieldSensitivityNode2);

                real tPhiA = mIntersection->mGeometry->evaluate_field_value(mIntersection->mNodes(0)->mNodeIndex,
                                                                                mIntersection->mNodes(0)->mCoordinates);
                real tPhiB = mIntersection->mGeometry->evaluate_field_value(mIntersection->mNodes(1)->mNodeIndex,
                                                                                mIntersection->mNodes(1)->mCoordinates);

                // Compute $\frac{\partial x_{\Gamma}}{\partial \phi}$
                Matrix<DDRMat> tDxgammaDphiA = -tPhiB / std::pow((tPhiA-tPhiB), 2)
                        * (mIntersection->mNodes(1)->mCoordinates - mIntersection->mNodes(0)->mCoordinates);
                Matrix<DDRMat> tDxgammaDphiB =  tPhiA / std::pow((tPhiA-tPhiB), 2)
                        * (mIntersection->mNodes(1)->mCoordinates - mIntersection->mNodes(0)->mCoordinates);

                // Compute dx/dp
                aSensitivities = trans(tDxgammaDphiA) * tFieldSensitivityNode1 + trans(tDxgammaDphiB) * tFieldSensitivityNode2;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
