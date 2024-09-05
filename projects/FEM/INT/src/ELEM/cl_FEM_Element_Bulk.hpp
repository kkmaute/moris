/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp"    //FEM/INT/src

namespace moris::fem
{
    class Set;
    //------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Element_Bulk : public Element
    {

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        Element_Bulk(){};

        /**
         * constructor
         * @param[ in ] aCell               a mesh cell pointer
         * @param[ in ] aSet                a fem set pointer
         * @param[ in ] aCluster            a fem cluster pointer
         * @param[ in ] aCellIndexInCluster an index for cell in cluster
         */
        Element_Bulk(
                mtk::Cell const   *aCell,
                Set               *aSet,
                Cluster           *aCluster,
                moris::moris_index aCellIndexInCluster );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_Bulk() override;

        //------------------------------------------------------------------------------
        /**
         * compute jacobian
         */
        void compute_jacobian() override;

        //------------------------------------------------------------------------------
        /**
         * compute residual
         */
        void compute_residual() override;

        //------------------------------------------------------------------------------
        /**
         * compute jacobian and residual
         */
        void compute_jacobian_and_residual() override;

        //------------------------------------------------------------------------------
        /**
         * compute dRdp
         */
        void compute_dRdp() override;

        //------------------------------------------------------------------------------
        /**
         * compute quantities of interest in a global way
         */
        void compute_QI() override;

        //------------------------------------------------------------------------------
        /**
         * compute dQIdp
         */
        void compute_dQIdp_explicit() override;

        //------------------------------------------------------------------------------
        /**
         * compute dRdp and dQIdp
         */
        void compute_dRdp_and_dQIdp() override;

        //------------------------------------------------------------------------------
        /**
         * compute dQIdu
         */
        void compute_dQIdu() override;

        //------------------------------------------------------------------------------
        /**
         * compute quantity of interest in a global way
         * @param[ in ] aFemMeshIndex mesh index to define on which mesh
         *                              values are evaluated
         */
        void compute_quantity_of_interest_global( const uint aFemMeshIndex ) override;

        //------------------------------------------------------------------------------
        /**
         * compute quantity of interest in a elemental way
         * @param[ in ] aFemMeshIndex mesh index to define on which mesh
         * @param[ in ] aAverageOutput flag, turn on to request the averaged, rather than integrated quantity on the element/facet
         */
        void compute_quantity_of_interest_elemental(
                const uint aFemMeshIndex,
                const bool aAverageOutput ) override;

        //------------------------------------------------------------------------------
        /**
         * compute quantity of interest in a elemental way
         * @param[ in ] Values
         */
        void compute_quantity_of_interest_elemental(
                Matrix< DDRMat > &tValues,
                uint              aIQIIndex,
                real             &aSpaceTimeVolume ) override;

        //------------------------------------------------------------------------------
        /**
         * compute volume over the element
         */
        real compute_volume( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) override;

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        /**
         * initialize the geometry interpolator for the IG element
         * using the mesh and GE
         * @param[ in ] aGeoLocalAssembly matrix with pdv local assembly indices
         *                               ( tNumNodes x tNumPdvTypes )
         */
        void init_ig_geometry_interpolator( Matrix< DDSMat > &aGeoLocalAssembly );

        //------------------------------------------------------------------------------
        /**
         * initialize the geometry interpolator for the IG element
         * using the mesh only
         */
        void init_ig_geometry_interpolator();
    };

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_ */
