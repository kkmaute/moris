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

        Intersection_Node::Intersection_Node(
                uint                      aFirstNodeIndex,
                uint                      aSecondNodeIndex,
                const Matrix<DDRMat>&     aFirstNodeCoordinates,
                const Matrix<DDRMat>&     aSecondNodeCoordinates,
                std::shared_ptr<Geometry> aInterfaceGeometry,
                real                      aIsocontourThreshold)
                : Child_Node({{aFirstNodeIndex, aSecondNodeIndex}},
                             {aFirstNodeCoordinates, aSecondNodeCoordinates},
                             xtk::Linear_Basis_Function(),
                             Interpolation::linear_interpolation_value(
                                     Matrix<DDRMat>(
                                             {{aInterfaceGeometry->get_field_value(aFirstNodeIndex, aFirstNodeCoordinates)},
                                              {aInterfaceGeometry->get_field_value(aSecondNodeIndex, aSecondNodeCoordinates)}}),
                                     aIsocontourThreshold)),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mFirstParentOnInterface(mInterfaceGeometry->get_field_value(aFirstNodeIndex, aFirstNodeCoordinates) == 0.0),
                  mSecondParentOnInterface(mInterfaceGeometry->get_field_value(aSecondNodeIndex, aSecondNodeCoordinates) == 0.0),
                  mGlobalCoordinates((mBasisValues(0) * aFirstNodeCoordinates) + (mBasisValues(1) * aSecondNodeCoordinates))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Intersection_Node::first_parent_on_interface()
        {
            return mFirstParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Intersection_Node::second_parent_on_interface()
        {
            return mSecondParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::set_starting_pdv_id(uint aStartingPdvIndex)
        {
            mStartingPdvIndex = aStartingPdvIndex;
            mPdvIndexSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Intersection_Node::get_starting_pdv_id()
        {
            //
            MORIS_ASSERT(mPdvIndexSet, "Starting PDV index must be set for an intersection.");
            return mStartingPdvIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Intersection_Node::get_coordinate_value(uint aCoordinateIndex)
        {
            return mGlobalCoordinates(aCoordinateIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_global_coordinates()
        {
            return mGlobalCoordinates;
        }

        moris::uint Intersection_Node::get_num_pdvs()
        {
            return mGlobalCoordinates.numel();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::get_all_sensitivities(Matrix<DDRMat>& aSensitivities)
        {
            // Geometry field sensitivities
            Matrix<DDRMat> tFieldSensitivityNode1(0, 0);
            Matrix<DDRMat> tFieldSensitivityNode2(0, 0);
            mInterfaceGeometry->get_field_adv_sensitivities(mParentNodeIndices(0),
                                                     mParentNodeCoordinates(0),
                                                     tFieldSensitivityNode1);
            mInterfaceGeometry->get_field_adv_sensitivities(mParentNodeIndices(1),
                                                     mParentNodeCoordinates(1),
                                                     tFieldSensitivityNode2);

            // Field values
            real tPhiA = mInterfaceGeometry->get_field_value(mParentNodeIndices(0),
                                                                  mParentNodeCoordinates(0));
            real tPhiB = mInterfaceGeometry->get_field_value(mParentNodeIndices(1),
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
