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
#include "cl_FEM_Element.hpp"               //FEM/INT/src

namespace moris
{
    namespace fem
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
                        mtk::Cell    const * aCell,
                        Set                * aSet,
                        Cluster            * aCluster,
                        moris::moris_index   aCellIndexInCluster );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~Element_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * compute jacobian
                 */
                void compute_jacobian();

                //------------------------------------------------------------------------------
                /**
                 * compute residual
                 */
                void compute_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute jacobian and residual
                 */
                void compute_jacobian_and_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp
                 */
                void compute_dRdp();

                //------------------------------------------------------------------------------
                /**
                 * compute quantities of interest in a global way
                 */
                void compute_QI();

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
                 * compute dQIdu
                 */
                void compute_dQIdu();

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a global way
                 * @param[ in ] aMeshIndex mesh index to define on which mesh
                 *                              values are evaluated
                 */
                void compute_quantity_of_interest_global( const uint aMeshIndex );

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a elemental way
                 * @param[ in ] aMeshIndex mesh index to define on which mesh
                 *                              values are evaluated
                 */
                void compute_quantity_of_interest_elemental( const uint aMeshIndex );

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a elemental way
                 * @param[ in ] Values
                 */
                void compute_quantity_of_interest_elemental(
                        Matrix< DDRMat > & tValues,
                        uint               aIQIIndex,
                        real             & aSpaceTimeVolume );

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
                 * using the mesh and GE
                 * @param[ in ] aGeoLocalAssembly matrix with pdv local assembly indices
                 *                               ( tNumNodes x tNumPdvTypes )
                 */
                void init_ig_geometry_interpolator( Matrix< DDSMat > & aGeoLocalAssembly );

                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG element
                 * using the mesh only
                 */
                void init_ig_geometry_interpolator();
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_ */

