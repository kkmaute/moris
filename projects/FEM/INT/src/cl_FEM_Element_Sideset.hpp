/*
 * cl_FEM_Element_Sideset.hpp
 *
 *  Created on: Mar 07, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

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
                        mtk::Cell const  * aCell,
                        Set              * aElementBlock,
                        Cluster          * aCluster,
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
                 * compute dRdp by analytical formulation
                 */
                void compute_dRdp();

                /**
                 * compute dRdp by finite difference
                 */
                void compute_dRdp_FD();

                //------------------------------------------------------------------------------
                /**
                 * compute quantities of interest
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
                        const std::string & aQINamee );

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
                 * @param[ in ] aSideOrdinal side ordinal for the master element
                 * @param[ in ] aIsActiveDv  list of if design variable is active
                 *                           (vertexIndex)(DvType) for master element
                 */
                void init_ig_geometry_interpolator(
                        uint                              aSideOrdinal,
                        moris::Cell< Matrix< DDSMat > > & aIsActiveDv );

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_ */
