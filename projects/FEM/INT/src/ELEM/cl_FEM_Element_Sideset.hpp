/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Sideset.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        class Set;

        //------------------------------------------------------------------------------

        /**
         * \brief Element_Sideset class
         */
        class Element_Sideset : public Element
        {
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * constructor
             * @param[ in ]     pointer to mesh interface object
             * @param[ in ]     cell of pointers to integrand of weak form of governing eqs.
             * @param[ in ]     cell of pointer to fem nodes
             * @param[ in ]     Pointer to element block
             */
            Element_Sideset(
                    mtk::Cell const   *aCell,
                    Set               *aElementBlock,
                    Cluster           *aCluster,
                    moris::moris_index aCellIndexInCluster );

            //------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Element_Sideset();

            //------------------------------------------------------------------------------

            /**
             * compute jacobian over the element
             */
            void compute_jacobian();

            //------------------------------------------------------------------------------
            /**
             * compute residual over the element
             */
            void compute_residual();

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
             * compute quantities of interest
             */
            void compute_QI();

            //------------------------------------------------------------------------------

            /**
             * compute dQIdu
             */
            void compute_dQIdu();

            //------------------------------------------------------------------------------

            /**
             * compute dQIdp
             */
            void compute_dQIdp_explicit();

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
             */
            real compute_volume( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------

          protected:
            //------------------------------------------------------------------------------

            /**
             * initialize the geometry interpolator for the IG element
             * @param[ in ] aSideOrdinal side ordinal for the leader element
             */
            void init_ig_geometry_interpolator( uint aSideOrdinal );

            //------------------------------------------------------------------------------

            /**
             * initialize the geometry interpolator for the IG element
             * @param[ in ] aSideOrdinal      side ordinal for the leader element
             * @param[ in ] aGeoLocalAssembly matrix with pdv local assembly indices
             *                                ( NumVertexIndices x NumPdvTypes )
             */
            void init_ig_geometry_interpolator(
                    uint              aSideOrdinal,
                    Matrix< DDSMat > &aGeoLocalAssembly );

            //------------------------------------------------------------------------------

        };    // class fem::Element_Sideset

        //------------------------------------------------------------------------------

    }    // namespace fem
}    // namespace moris

#endif /* SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_ */
