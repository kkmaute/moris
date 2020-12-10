#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(
                uint                       aFirstParentNodeIndex,
                uint                       aSecondParentNodeIndex,
                const Matrix<DDRMat>&      aFirstParentNodeCoordinates,
                const Matrix<DDRMat>&      aSecondParentNodeCoordinates,
                Matrix<DDUMat>             aAncestorNodeIndices,
                Cell<Matrix<DDRMat>>       aAncestorNodeCoordinates,
                const xtk::Basis_Function& aAncestorBasisFunction,
                Matrix<DDRMat>             aAncestorLocalCoordinates,
                std::shared_ptr<Geometry>  aInterfaceGeometry,
                real                       aIsocontourThreshold,
                real                       aIsocontourTolerance)
                : Child_Node(
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        aAncestorLocalCoordinates),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mFirstParentOnInterface( std::abs(
                          mInterfaceGeometry->get_field_value(aFirstParentNodeIndex, aFirstParentNodeCoordinates) -
                          aIsocontourThreshold) < aIsocontourTolerance ),
                  mSecondParentOnInterface( std::abs(
                          mInterfaceGeometry->get_field_value(aSecondParentNodeIndex, aSecondParentNodeCoordinates) -
                          aIsocontourThreshold) < aIsocontourTolerance ),
                  mGlobalCoordinates(
                          mBasisValues(0) * aFirstParentNodeCoordinates + mBasisValues(1) * aSecondParentNodeCoordinates)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Intersection_Node::get_ancestor_coordinate_determining_adv_ids(uint aAncestorIndex)
        {
            return mInterfaceGeometry->get_determining_adv_ids(
                    mAncestorNodeIndices(aAncestorIndex), mAncestorNodeCoordinates(aAncestorIndex));
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

        real Intersection_Node::get_local_coordinate()
        {
            return mLocalCoordinate;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_global_coordinates()
        {
            return mGlobalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Intersection_Node::get_coordinate_value(uint aCoordinateIndex)
        {
            return mGlobalCoordinates(aCoordinateIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Intersection_Node::get_num_pdvs()
        {
            return mGlobalCoordinates.numel();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::set_starting_pdv_id(moris_id aPDVStartingID)
        {
            mPDVStartingID = aPDVStartingID;
            mPDVStartingIDSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Intersection_Node::get_starting_pdv_id()
        {
            MORIS_ASSERT(mPDVStartingIDSet, "PDV Starting ID must be set for an intersection.");
            return mPDVStartingID;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::set_id(moris_id aNodeID)
        {
            mNodeID = aNodeID;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Intersection_Node::set_owner(moris_index aNodeOwner)
        {
            mNodeOwner = aNodeOwner;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Intersection_Node::get_id()
        {
            return mNodeID;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index Intersection_Node::get_owner()
        {
            return mNodeOwner;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
