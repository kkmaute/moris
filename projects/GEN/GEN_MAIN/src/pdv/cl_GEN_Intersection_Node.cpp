#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(
                real                       aLocalCoordinate,
                real                       aFirstParentNodeIndex,
                real                       aSecondParentNodeIndex,
                const Matrix<DDRMat>&      aFirstParentNodeLocalCoordinates,
                const Matrix<DDRMat>&      aSecondParentNodeLocalCoordinates,
                Matrix<DDUMat>             aAncestorNodeIndices,
                Cell<Matrix<DDRMat>>       aAncestorNodeCoordinates,
                const xtk::Basis_Function& aAncestorBasisFunction,
                std::shared_ptr<Geometry>  aInterfaceGeometry,
                real                       aIsocontourThreshold,
                real                       aIsocontourTolerance,
                real                       aIntersectionTolerance)
                : Child_Node(
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        0.5 * ((1 - aLocalCoordinate) * aFirstParentNodeLocalCoordinates + (1 + aLocalCoordinate) * aSecondParentNodeLocalCoordinates)),
                  mLocalCoordinate(aLocalCoordinate),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mIsocontourThreshold(aIsocontourThreshold)
        {
            // Parent basis
            Matrix<DDRMat> tFirstParentBasisValues;
            Matrix<DDRMat> tSecondParentBasisValues;
            aAncestorBasisFunction.evaluate_basis_function(aFirstParentNodeLocalCoordinates, tFirstParentBasisValues);
            aAncestorBasisFunction.evaluate_basis_function(aSecondParentNodeLocalCoordinates, tSecondParentBasisValues);

            // Global coordinates of intersection and parents
            mGlobalCoordinates = mBasisValues(0) * aAncestorNodeCoordinates(0);
            Matrix<DDRMat> tFirstParentGlobalCoordinates = tFirstParentBasisValues(0) * aAncestorNodeCoordinates(0);
            Matrix<DDRMat> tSecondParentGlobalCoordinates = tSecondParentBasisValues(0) * aAncestorNodeCoordinates(0);
            for (uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++)
            {
                mGlobalCoordinates += mBasisValues(tBasisIndex) * aAncestorNodeCoordinates(tBasisIndex);
                tFirstParentGlobalCoordinates += tFirstParentBasisValues(tBasisIndex) * aAncestorNodeCoordinates(tBasisIndex);
                tSecondParentGlobalCoordinates += tSecondParentBasisValues(tBasisIndex) * aAncestorNodeCoordinates(tBasisIndex);
            }

            // Parents on interface
            real tFirstParentPhi = aInterfaceGeometry->get_field_value(aFirstParentNodeIndex, tFirstParentGlobalCoordinates);
            real tSecondParentPhi = aInterfaceGeometry->get_field_value(aSecondParentNodeIndex, tSecondParentGlobalCoordinates);
            real tParentLength = norm(tSecondParentGlobalCoordinates - tFirstParentGlobalCoordinates);
            mFirstParentOnInterface = std::abs(tFirstParentPhi) < aIsocontourTolerance or
                    0.5 * tParentLength * std::abs(1 + aLocalCoordinate) < aIntersectionTolerance;
            mSecondParentOnInterface = std::abs(tSecondParentPhi) < aIsocontourTolerance or
                    0.5 * tParentLength * std::abs(1 - aLocalCoordinate) < aIntersectionTolerance;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_dcoordinate_dadv_from_ancestor(uint aAncestorIndex)
        {
            MORIS_ASSERT(aAncestorIndex < this->get_num_ancestors(),
                    "Ancestor index for intersection sensitivities out of bounds for the number of ancestors.");

            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get geometry field sensitivity with respect to ADVs
            const Matrix<DDRMat>& tFieldSensitivities = tLockedInterfaceGeometry->get_field_sensitivities(
                    mAncestorNodeIndices(aAncestorIndex),
                    mAncestorNodeCoordinates(aAncestorIndex));

            // Compute full sensitivity of global coordinates with respect to ADVs
            return 0.5 * this->get_dcoordinate_dfield_from_ancestor(aAncestorIndex) *
                    trans(mAncestorNodeCoordinates(1) - mAncestorNodeCoordinates(0)) * tFieldSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Intersection_Node::get_ancestor_coordinate_determining_adv_ids(uint aAncestorIndex)
        {
            MORIS_ASSERT(aAncestorIndex < this->get_num_ancestors(),
                    "Ancestor index for intersection determining ADV IDs out of bounds for the number of ancestors.");

            return mInterfaceGeometry.lock()->get_determining_adv_ids(
                    mAncestorNodeIndices(aAncestorIndex), mAncestorNodeCoordinates(aAncestorIndex));
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Intersection_Node::parent_edge_is_intersected()
        {
            return (std::abs(mLocalCoordinate) <= 1.0);
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
