/*
 * cl_FEM_Element_Bulk.hpp
 *
 *  Created on: Apr 22, 2018
 *      Author: schmidt
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
                 * compute dRdp by analytical formulation
                 */
                void compute_dRdp();

                /**
                 * compute dRdp by finite difference
                 */
                void compute_dRdp_FD();

                //------------------------------------------------------------------------------
                /**
                 * compute quantities of interest in a global way
                 */
                void compute_QI();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp by analytical formulation
                 */
                void compute_dQIdp_explicit();

                /**
                 * compute dQIdp by finite difference
                 */
                void compute_dQIdp_explicit_FD();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdu
                 */
                void compute_dQIdu();

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a global way
                 * @param[ in ] aOutputType an enum for the output type
                 */
                void compute_quantity_of_interest_global(
                        const uint          aMeshIndex,
                        const std::string & aQIName );

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a nodal way
                 * @param[ in ] aOutputType an enum for the output type
                 */
                void compute_quantity_of_interest_nodal(
                        const uint          aMeshIndex,
                        const std::string & aQIName );

                //------------------------------------------------------------------------------
                /**
                 * compute quantity of interest in a elemental way
                 * @param[ in ] aOutputType an enum for the output type
                 */
                void compute_quantity_of_interest_elemental(
                        const uint          aMeshIndex,
                        const std::string & aQIName );

                //------------------------------------------------------------------------------
                /**
                 * compute volume over the element
                 */
                real compute_volume( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
            protected:
                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG element
                 * @param[ in ] aIsActiveDv list of if design variable is active
                 *                          (vertexIndex)(DvType) for master element
                 */
                void init_ig_geometry_interpolator(
                        moris::Cell< Matrix< DDSMat > > & aIsActiveDv );
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_ */
