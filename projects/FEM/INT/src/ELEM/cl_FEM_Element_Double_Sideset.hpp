/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Double_Sideset.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp"    //FEM/INT/src

namespace moris::fem
{
    class Set;

    //------------------------------------------------------------------------------
    /**
     * @brief Element_Double_Sideset class
     */
    class Element_Double_Sideset : public Element
    {

      public:
        //------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aLeaderIGCell         pointer to mesh cell for leader element
         * @param[ in ] aFollowerIGCell        pointer to mesh cell for follower element
         * @param[ in ] aSet                pointer to FEM set to which the elements belong
         * @param[ in ] aCluster            pointer to FEM cluster to which the elements belong
         * @param[ in ] aCellIndexInCluster index of the element in the cluster
         */
        Element_Double_Sideset(
                mtk::Cell const         *aLeaderIGCell,
                mtk::Cell const         *aFollowerIGCell,
                Set                     *aSet,
                Cluster                 *aCluster,
                moris::moris_index const aCellIndexInCluster );

        //------------------------------------------------------------------------------

        ~Element_Double_Sideset() override = default;

        //------------------------------------------------------------------------------
        /**
         * compute residual over the element
         */
        void compute_residual() override;

        //------------------------------------------------------------------------------
        /**
         * compute jacobian over the element
         */
        void compute_jacobian() override;

        /**
         * compute jacobian and residual over the element
         */
        void compute_jacobian_and_residual() override;

        //------------------------------------------------------------------------------

        /**
         * compute dRdp
         */
        void compute_dRdp() override;

        //------------------------------------------------------------------------------

        /**
         * compute QI
         */
        void compute_QI() override;

        //------------------------------------------------------------------------------

        /**
         * compute dQIdu
         */
        void
        compute_dQIdu() override
        {
            MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdu - not implemented." );
        }

        //------------------------------------------------------------------------------

        /**
         * compute dQIdp
         */
        void
        compute_dQIdp_explicit() override
        {
            MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdp_explicit - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute dRdp and dQIdp
         */
        void compute_dRdp_and_dQIdp() override;

        //------------------------------------------------------------------------------

        /**
         * compute quantity of interest in a global way
         * @param[ in ] aFemMeshIndex mesh index for used IG mesh
         */
        void compute_quantity_of_interest_global( const uint aFemMeshIndex ) override;

        //------------------------------------------------------------------------------

        /**
         * compute quantity of interest in a elemental way
         * @param[ in ] aFemMeshIndex mesh index for used IG mesh
         * @param[ in ] aAverageOutput flag, turn on to request the averaged, rather than integrated quantity on the element/facet
         */
        void compute_quantity_of_interest_elemental(
                const uint aFemMeshIndex,
                const bool aAverageOutput ) override;

        //------------------------------------------------------------------------------

        /**
         * compute volume over the element
         * @param[ in ] aIsLeader enum for leader or follower
         */
        real compute_volume( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) override;

        //------------------------------------------------------------------------------

        void init_ig_geometry_interpolator() const;

        //------------------------------------------------------------------------------

      protected:
        /**
         * initialize the geometry interpolator for the IG leader and follower element
         * @param[ in ] aLeaderSideOrdinal side ordinal for the leader element
         * @param[ in ] aFollowerSideOrdinal  side ordinal for the follower element
         */
        virtual void initialize_leader_follower_ig_interpolator( mtk::Leader_Follower const aLeaderFollowerType ) const;

        //------------------------------------------------------------------------------
        Matrix< DDSMat > get_local_cluster_assembly_indices(
                moris_index const aLeaderSideOrdinal,
                moris_index const aFollowerSideOrdinal ) const;

        //------------------------------------------------------------------------------
        /**
         * @brief Return the index at which the leader cell is stored in the vector of primary cells in the parent cluster.
         * @details If the cluster contains cells [1, 3, 4, 7] and the cell is 4, then the index is 2.
         * @return Index of the cell in the vector of primary cells in the parent cluster.
         */
        moris_index get_leader_local_cell_index() const;

        //------------------------------------------------------------------------------
        /**
         * @brief Return the index at which the follower cell is stored in the vector of primary cells in the parent cluster.
         * @attention For the double sideset element, this is the same as the leader cell index.
         * This is not the case for the nonconformal sideset element which is derived from this class.
         * @details If the cluster contains cells [1, 3, 4, 7] and the cell is 4, then the index is 2.
         * @return Index of the cell in the vector of primary cells in the parent cluster.
         */
        virtual moris_index get_follower_local_cell_index() const;

        //------------------------------------------------------------------------------

        virtual Matrix< DDRMat > get_leader_normal( uint aGPIndex ) const
        {
            moris_index const tLeaderSideOrd = mCluster->mLeaderListOfSideOrdinals( get_leader_local_cell_index() );
            return mCluster->get_side_normal( mLeaderCell, tLeaderSideOrd );
        }

      private:
        Matrix< DDRMat > get_follower_integration_point( uint aGPIndex ) const override;

        //------------------------------------------------------------------------------
        /**
         * @brief Get the leader and follower quadrature point coordinates and set the field interpolators to this coordinate
         */
        void initialize_quadrature_point( uint iGP );
    };

}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_ */
