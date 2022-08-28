/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Time_Boundary.hpp
 *
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
                void compute_dQIdu();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                void compute_dQIdp_explicit();

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
                 */
                void init_ig_geometry_interpolator( uint aTimeOrdinal );

                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG element
                 * @param[ in ] aTimeOrdinal      time ordinal for the IG element
                 * @param[ in ] aGeoLocalAssembly matrix to fill with pdv local assembly indices
                 *                               ( NumNodes x NumPdvTypes )
                 */
                void init_ig_geometry_interpolator(
                        uint               aTimeOrdinal,
                        Matrix< DDSMat > & aGeoLocalAssembly );

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_Element_Time_Boundary_HPP_ */

