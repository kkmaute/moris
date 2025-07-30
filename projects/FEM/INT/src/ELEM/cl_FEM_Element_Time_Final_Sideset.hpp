/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Time_Final_Sideset.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_Element_Time_Final_Sideset_HPP_
#define SRC_FEM_CL_FEM_Element_Time_Final_Sideset_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris::fem
{
    class Set;
    //------------------------------------------------------------------------------
    /**
     * \brief Element_Time_Final_Sideset class
     */
    class Element_Time_Final_Sideset : public Element
    {

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        real mTimeTol   = 1.0e-6;
        bool mSkipTime  = false;

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
        Element_Time_Final_Sideset(
                mtk::Cell const   *aCell,
                Set               *aSet,
                Cluster           *aCluster,
                moris::moris_index aCellIndexInCluster );

        //------------------------------------------------------------------------------
        /**
         * destructor
         */
            ~Element_Time_Final_Sideset() override;

            //------------------------------------------------------------------------------

            /**
             * set final time
             * @param[ in ]     final time
             */
            void set_final_time( real aFinalTime );

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

            //------------------------------------------------------------------------------
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
            void compute_dQIdu() override;

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp
             */
            void compute_dQIdp_explicit() override;

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp by finite difference
             */
            void compute_dRdp_and_dQIdp() override;

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * initialize the geometry interpolator for the IG element
             */
            void init_ig_geometry_interpolator();

            //------------------------------------------------------------------------------
            /**
             * initialize the geometry interpolator for the IG element
             * @param[ in ] aGeoLocalAssembly matrix with pdv local assembly indices
             *                                ( tNumNodes x tNumPdvTypes )
             */
            void init_ig_geometry_interpolator( Matrix< DDSMat > &aGeoLocalAssembly );

            //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_Element_Time_Final_Sideset_HPP_ */
