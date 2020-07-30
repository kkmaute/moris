/*
 * cl_FEM_Element_Time_Boundary.hpp
 *
 *  Created on: Mar 19, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_Element_Time_Boundary_HPP_
#define SRC_FEM_CL_FEM_Element_Time_Boundary_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
        class Set;
        //------------------------------------------------------------------------------
        /**
         * \brief Element_Time_Boundary class
         */
        class Element_Time_Boundary : public Element
        {

                //------------------------------------------------------------------------------
            protected:
                //------------------------------------------------------------------------------


                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------

                /**
                 * constructor
                 *
                 * @param[ in ]     pointer to mesh interface object
                 * @param[ in ]     cell of pointers to integrand of weak form of governing eqs.
                 * @param[ in ]     cell of pointer to fem nodes
                 * @param[ in ]     Pointer to element block
                 */
                Element_Time_Boundary(
                        mtk::Cell const  * aCell,
                        Set              * aSet,
                        Cluster          * aCluster,
                        moris::moris_index aCellIndexInCluster );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Element_Time_Boundary();

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
                 * compute volume over the element
                 */
                real compute_volume( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
                {
                    MORIS_ERROR( false, "Element_Time_Boundary::compute_volume - not implemented." );
                    return 0.0;
                }

                //------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG element
                 * @param[ in ] aTimeOrdinal time ordinal for the IG element
                 * @param[ in ] aIsActiveDv  list of if design variable is active
                 *                           (vertexIndex)(DvType) for master element
                 */
                void init_ig_geometry_interpolator(
                        uint                              aTimeOrdinal,
                        moris::Cell< Matrix< DDSMat > > & aIsActiveDv );

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_Element_Time_Boundary_HPP_ */
