#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_XTK_Linear_Basis_Functions.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(Matrix<DDUMat>            aParentNodeIndices,
                                             Cell<Matrix<DDRMat>>      aParentNodeCoordinates,
                                             std::shared_ptr<Geometry> aInterfaceGeometry,
                                             real                      aIsocontourThreshold)
                : Child_Node(aParentNodeIndices,
                             aParentNodeCoordinates,
                             xtk::Linear_Basis_Function(),
                             Interpolation::linear_interpolation_value(
                                     Matrix<DDRMat>(
                                             {{aInterfaceGeometry->evaluate_field_value(aParentNodeIndices(0), aParentNodeCoordinates(0))},
                                              {aInterfaceGeometry->evaluate_field_value(aParentNodeIndices(1), aParentNodeCoordinates(1))}}),
                                     aIsocontourThreshold)),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mCoordinates((mBasisValues(0) * aParentNodeCoordinates(0)) + (mBasisValues(1) * aParentNodeCoordinates(1)))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::set_starting_pdv_index(uint aStartingPdvIndex)
        {
            mStartingPdvIndex = aStartingPdvIndex;
            mPdvIndexSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Intersection_Node::get_starting_pdv_index()
        {
            MORIS_ASSERT(mPdvIndexSet, "Starting PDV index must be set for an intersection.");
            return mStartingPdvIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Intersection_Node::get_coordinate_value(uint aCoordinateIndex)
        {
            return mCoordinates(aCoordinateIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::get_all_sensitivities(Matrix<DDRMat>& aSensitivities)
        {
            MORIS_ERROR(mInterfaceGeometry->sensitivities_available(),
                        "Sensitivities of an intersection not implemented if a geometry does not have its sensitivities available.");

            // Geometry field sensitivities
            Matrix<DDRMat> tFieldSensitivityNode1(0, 0);
            Matrix<DDRMat> tFieldSensitivityNode2(0, 0);
            mInterfaceGeometry->evaluate_sensitivity(mParentNodeIndices(0),
                                                     mParentNodeCoordinates(0),
                                                     tFieldSensitivityNode1);
            mInterfaceGeometry->evaluate_sensitivity(mParentNodeIndices(1),
                                                     mParentNodeCoordinates(1),
                                                     tFieldSensitivityNode2);

            // Field values
            real tPhiA = mInterfaceGeometry->evaluate_field_value(mParentNodeIndices(0),
                                                                  mParentNodeCoordinates(0));
            real tPhiB = mInterfaceGeometry->evaluate_field_value(mParentNodeIndices(1),
                                                                  mParentNodeCoordinates(1));

            // Compute $\frac{\partial x_{\Gamma}}{\partial \phi}$
            Matrix<DDRMat> tDxgammaDphiA = -tPhiB / std::pow((tPhiA - tPhiB), 2)
                    * (mParentNodeCoordinates(1) - mParentNodeCoordinates(0));
            Matrix<DDRMat> tDxgammaDphiB =  tPhiA / std::pow((tPhiA - tPhiB), 2)
                    * (mParentNodeCoordinates(1) - mParentNodeCoordinates(0));

            // Compute dx/dp
            aSensitivities = trans(tDxgammaDphiA) * tFieldSensitivityNode1 + trans(tDxgammaDphiB) * tFieldSensitivityNode2;
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}