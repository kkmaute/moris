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

namespace moris
{
    namespace fem
    {
        class Set;
        //------------------------------------------------------------------------------
        /**
         * @brief Element_Double_Sideset class
         */
        class Element_Double_Sideset : public Element
        {

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * constructor
             * @param[ in ] aLeftIGCell         pointer to mesh cell for leader element
             * @param[ in ] aRightIGCell        pointer to mesh cell for follower element
             * @param[ in ] aSet                pointer to FEM set to which the elements belong
             * @param[ in ] aCluster            pointer to FEM cluster to which the elements belong
             * @param[ in ] aCellIndexInCluster index of the element in the cluster
             */
            Element_Double_Sideset(
                    mtk::Cell const   *aLeftIGCell,
                    mtk::Cell const   *aRightIGCell,
                    Set               *aSet,
                    Cluster           *aCluster,
                    moris::moris_index aCellIndexInCluster );

            //------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Element_Double_Sideset();

            //------------------------------------------------------------------------------

            /**
             * compute residual over the element
             */
            void compute_residual();

            //------------------------------------------------------------------------------

            /**
             * compute jacobian over the element
             */
            void compute_jacobian();

            //------------------------------------------------------------------------------

            /**
             * compute jacobian and residual over the element
             */
            void compute_jacobian_and_residual();

            //------------------------------------------------------------------------------

            /**
             * compute dRdp
             */
            void compute_dRdp();

            //------------------------------------------------------------------------------

            /**
             * compute QI
             */
            void compute_QI();

            //------------------------------------------------------------------------------

            /**
             * compute dQIdu
             */
            void
            compute_dQIdu()
            {
                MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdu - not implemented." );
            }

            //------------------------------------------------------------------------------

            /**
             * compute dQIdp
             */
            void
            compute_dQIdp_explicit()
            {
                MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdp_explicit - not implemented." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dRdp and dQIdp
             */
            void compute_dRdp_and_dQIdp();

            //------------------------------------------------------------------------------

            /**
             * compute quantity of interest in a global way
             * @param[ in ] aFemMeshIndex mesh index for used IG mesh
             */
            void compute_quantity_of_interest_global( const uint aFemMeshIndex );

            //------------------------------------------------------------------------------

            /**
             * compute quantity of interest in a elemental way
             * @param[ in ] aFemMeshIndex mesh index for used IG mesh
             * @param[ in ] aAverageOutput flag, turn on to request the averaged, rather than integrated quantity on the element/facet
             */
            void compute_quantity_of_interest_elemental(
                    const uint aFemMeshIndex,
                    const bool aAverageOutput );

            //------------------------------------------------------------------------------

            /**
             * compute volume over the element
             * @param[ in ] aIsLeader enum for leader or follower
             */
            real compute_volume( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------

          protected:
            //------------------------------------------------------------------------------

            /**
             * initialize the geometry interpolator for the IG leader and follower element
             * @param[ in ] aLeaderSideOrdinal side ordinal for the leader element
             * @param[ in ] aFollowerSideOrdinal  side ordinal for the follower element
             */
            void init_ig_geometry_interpolator(
                    uint aLeaderSideOrdinal,
                    uint aFollowerSideOrdinal );

            /**
             * @brief Return the index at which the leader cell is stored in the vector of primary cells in the parent cluster.
             * @details If the cluster contains cells [1, 3, 4, 7] and the cell is 4, then the index is 2.
             * @return Index of the cell in the vector of primary cells in the parent cluster.
             */
            moris_index get_leader_local_cell_index() const;

            /**
             * @brief Return the index at which the follower cell is stored in the vector of primary cells in the parent cluster.
             * @attention For the double sideset element, this is the same as the leader cell index.
             * This is not the case for the nonconformal sideset element which is derived from this class.
             * @details If the cluster contains cells [1, 3, 4, 7] and the cell is 4, then the index is 2.
             * @return Index of the cell in the vector of primary cells in the parent cluster.
             */
            virtual moris_index get_follower_local_cell_index() const;

            //------------------------------------------------------------------------------

            /**
             * initialize the geometry interpolator for the IG leader and follower element
             * @param[ in ] aLeaderSideOrdinal side ordinal for the leader element
             * @param[ in ] aFollowerSideOrdinal  side ordinal for the follower element
             * @param[ in ] aGeoLocalAssembly  matrix with pdv local assembly indices
             *                                 for leader element
             *                                 ( NumVertexIndices x NumPdvTypes )
             */
            void init_ig_geometry_interpolator(
                    uint              aLeaderSideOrdinal,
                    uint              aFollowerSideOrdinal,
                    Matrix< DDSMat > &aGeoLocalAssembly );

          private:
            Matrix< DDRMat > get_follower_integration_point( uint aGPIndex ) const override;
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_ */
