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
                  mFirstParentOnInterface( std::abs( mInterfaceGeometry->get_field_value(aFirstNodeIndex, aFirstNodeCoordinates) ) <=  1E-10),
                  mSecondParentOnInterface(std::abs(mInterfaceGeometry->get_field_value(aSecondNodeIndex, aSecondNodeCoordinates) ) <= 1E-10),
                  mGlobalCoordinates((mBasisValues(0) * aFirstNodeCoordinates) + (mBasisValues(1) * aSecondNodeCoordinates))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(
                uint                      aFirstNodeIndex,
                uint                      aSecondNodeIndex,
                const Matrix<DDRMat>&     aFirstNodeCoordinates,
                const Matrix<DDRMat>&     aSecondNodeCoordinates,
                std::shared_ptr<Geometry> aInterfaceGeometry)
                : Child_Node({{aFirstNodeIndex, aSecondNodeIndex}},
                             {aFirstNodeCoordinates, aSecondNodeCoordinates},
                             xtk::Linear_Basis_Function(),
                             Matrix<DDRMat>( { {0.0} } ) ),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mFirstParentOnInterface(false),
                  mSecondParentOnInterface(false),
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

        uint Intersection_Node::get_num_pdvs()
        {
            return mGlobalCoordinates.numel();
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

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_first_parent_sensitivities()
        {
            // Get geometry field values
            real tPhi1 = mInterfaceGeometry->get_field_value(mParentNodeIndices(0), mParentNodeCoordinates(0));
            real tPhi2 = mInterfaceGeometry->get_field_value(mParentNodeIndices(1), mParentNodeCoordinates(1));

            // Get geometry field sensitivity with respect to ADVs
            const Matrix<DDRMat>& tFieldSensitivity = mInterfaceGeometry->get_field_sensitivities(
                    mParentNodeIndices(0),
                    mParentNodeCoordinates(0));

            // Compute sensitivity of the global coordinate with respect to the field value
            Matrix<DDRMat> tCoordinateSensitivity = -tPhi2 / std::pow((tPhi1 - tPhi2), 2)
                    * (mParentNodeCoordinates(1) - mParentNodeCoordinates(0));

            // Compute full sensitivity of global coordinates with respect to ADVs
            return (trans(tCoordinateSensitivity) * tFieldSensitivity);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_second_parent_sensitivities()
        {
            // Get geometry field values
            real tPhi1 = mInterfaceGeometry->get_field_value(mParentNodeIndices(0), mParentNodeCoordinates(0));
            real tPhi2 = mInterfaceGeometry->get_field_value(mParentNodeIndices(1), mParentNodeCoordinates(1));

            // Get geometry field sensitivity with respect to ADVs
            const Matrix<DDRMat>& tFieldSensitivity = mInterfaceGeometry->get_field_sensitivities(
                    mParentNodeIndices(1),
                    mParentNodeCoordinates(1));

            // Compute sensitivity of the global coordinate with respect to the field value
            Matrix<DDRMat> tCoordinateSensitivity = tPhi1 / std::pow((tPhi1 - tPhi2), 2)
                    * (mParentNodeCoordinates(1) - mParentNodeCoordinates(0));

            // Compute full sensitivity of global coordinates with respect to ADVs
            return (trans(tCoordinateSensitivity) * tFieldSensitivity);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Intersection_Node::get_first_parent_determining_adv_ids()
        {
            return mInterfaceGeometry->get_determining_adv_ids(mParentNodeIndices(0), mParentNodeCoordinates(0));
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Intersection_Node::get_second_parent_determining_adv_ids()
        {
            return mInterfaceGeometry->get_determining_adv_ids(mParentNodeIndices(1), mParentNodeCoordinates(1));
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}
