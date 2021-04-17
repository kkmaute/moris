#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(
                real                               aLocalCoordinate,
                std::shared_ptr<Intersection_Node> aFirstParentNode,
                std::shared_ptr<Intersection_Node> aSecondParentNode,
                real                               aFirstParentNodeIndex,
                real                               aSecondParentNodeIndex,
                const Matrix<DDRMat>&              aFirstParentNodeLocalCoordinates,
                const Matrix<DDRMat>&              aSecondParentNodeLocalCoordinates,
                Matrix<DDUMat>                     aAncestorNodeIndices,
                Cell<Matrix<DDRMat>>               aAncestorNodeCoordinates,
                const xtk::Basis_Function&         aAncestorBasisFunction,
                std::shared_ptr<Geometry>          aInterfaceGeometry,
                real                               aIsocontourThreshold,
                real                               aIsocontourTolerance,
                real                               aIntersectionTolerance)
                : Child_Node(
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        0.5 * ((1 - aLocalCoordinate) * aFirstParentNodeLocalCoordinates + (1 + aLocalCoordinate) * aSecondParentNodeLocalCoordinates)),
                  mLocalCoordinate(aLocalCoordinate),
                  mIsIntersected(false),
                  mInterfaceGeometry(aInterfaceGeometry),
                  mIsocontourThreshold(aIsocontourThreshold),
                  mFirstParentNode(aFirstParentNode),
                  mSecondParentNode(aSecondParentNode)
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
            

            if (mFirstParentOnInterface or mSecondParentOnInterface)
            {
                mIsIntersected = true;
            }
            else
            {
                mIsIntersected = (std::abs(mLocalCoordinate) <= 1.0);
            }

        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Intersection_Node::get_dcoordinate_dadv()
        {
            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get sensitivity values from other ancestors
            Matrix<DDRMat> tAncestorSensitivities;
            for (uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++)
            {
                // Get geometry field sensitivity with respect to ADVs
                const Matrix<DDRMat>& tFieldSensitivities = tLockedInterfaceGeometry->get_field_sensitivities(
                        mAncestorNodeIndices(tAncestorNode),
                        mAncestorNodeCoordinates(tAncestorNode));

                // Ancestor sensitivities
                tAncestorSensitivities = 0.5 * this->get_dcoordinate_dfield_from_ancestor(tAncestorNode) *
                        trans( trans(mAncestorNodeCoordinates(1) - mAncestorNodeCoordinates(0)) * tFieldSensitivities );

                // Join sensitivities
                uint tJoinedSensitivityLength = mCoordinateSensitivities.n_rows();
                mCoordinateSensitivities.resize(
                        tJoinedSensitivityLength + tAncestorSensitivities.n_rows(),
                        tAncestorSensitivities.n_cols());
                for (uint tAncestorSensitivity = 0; tAncestorSensitivity < tAncestorSensitivities.n_rows(); tAncestorSensitivity++)
                {
                    for (uint tCoordinateIndex = 0; tCoordinateIndex < tAncestorSensitivities.n_cols(); tCoordinateIndex++)
                    {
                        mCoordinateSensitivities(tJoinedSensitivityLength + tAncestorSensitivity, tCoordinateIndex) =
                                tAncestorSensitivities(tAncestorSensitivity, tCoordinateIndex);
                    }
                }
            }

            return mCoordinateSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Intersection_Node::get_coordinate_determining_adv_ids()
        {
            // Locked interface geometry
            std::shared_ptr<Geometry> tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get sensitivity values from other ancestors
            for (uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++)
            {
                // Get geometry field sensitivity with respect to ADVs
                const Matrix<DDSMat>& tAncestorADVIDs = tLockedInterfaceGeometry->get_determining_adv_ids(
                        mAncestorNodeIndices(tAncestorNode),
                        mAncestorNodeCoordinates(tAncestorNode));

                // Join sensitivities
                uint tJoinedSensitivityLength = mCoordinateDeterminingADVIDs.n_rows();
                mCoordinateDeterminingADVIDs.resize(tJoinedSensitivityLength + tAncestorADVIDs.n_rows(), 1);
                for (uint tAncestorSensitivity = 0; tAncestorSensitivity < tAncestorADVIDs.n_rows(); tAncestorSensitivity++)
                {
                    mCoordinateDeterminingADVIDs(tJoinedSensitivityLength + tAncestorSensitivity) =
                            tAncestorADVIDs(tAncestorSensitivity);
                }
            }

            return mCoordinateDeterminingADVIDs;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Intersection_Node::parent_edge_is_intersected()
        {
            return mIsIntersected;
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
