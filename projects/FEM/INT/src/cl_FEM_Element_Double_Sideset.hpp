/*
 * cl_FEM_Element_Double_Sideset.hpp
 *
 *  Created on: May 13, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_

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
        class Element_Double_Sideset : public Element
        {

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ]     pointer to mesh cell
                 * @param[ in ]     pointer to element block
                 * @param[ in ]     pointer to cluster
                 *
                 */
                Element_Double_Sideset(
                        mtk::Cell const  * aLeftIGCell,
                        mtk::Cell const  * aRightIGCell,
                        Set              * aSet,
                        Cluster          * aCluster,
                        moris::moris_index aCellIndexInCluster);

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Element_Double_Sideset();

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
                 * compute dRdp over the element
                 */
                void compute_dRdp();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                void compute_dQIdp_explicit();

                //------------------------------------------------------------------------------
                /**
                 * compute volume over the element
                 */
                real compute_volume( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
            protected:

                //------------------------------------------------------------------------------
                /**
                 * initialize integration geometry interpolator
                 */
                void init_ig_geometry_interpolator(
                        uint aMasterSideOrdinal,
                        uint aSlaveSideOrdinal );
                void init_ig_geometry_interpolator_with_pdv(
                        uint aMasterSideOrdinal,
                        uint aSlaveSideOrdinal,
                        moris::Cell< Matrix< DDSMat > > & aMasterIsActiveDv,
                        moris::Cell< Matrix< DDSMat > > & aSlaveIsActiveDv );
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_ */
