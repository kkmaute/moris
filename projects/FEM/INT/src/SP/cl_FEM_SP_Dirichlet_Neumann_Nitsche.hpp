/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Dirichlet_Neumann_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_DIRICHLET_NEUMANN_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_DIRICHLET_NEUMANN_NITSCHE_HPP_

#include <map>

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CON/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /*
         * Stabilization parameter for SlipBoundary BC on velocity for fluid problem
         * applied with Nitsche's formulation
         *
         * by Juntunen abd Stenberg 2009
         *
         * tangential direction:
         *             alpha_tang = 1/gamma ;
         *             gamma_t1 = alpha_tgang/(alpha_tang*sliplength+h)
         *             gamma_t2 = h/( alpha_tang*sliplength+h )
         */
        class SP_Dirichlet_Neumann_Nitsche : public Stabilization_Parameter
        {

            //------------------------------------------------------------------------------

          private:
            // populate the dof map (default)
            MSI::Dof_Type mMasterDofTemp = MSI::Dof_Type::THETA;

            // default tuple for element size to define cluster measure
            std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Master_Slave > mElementSizeTuple =
                    std::make_tuple(
                            fem::Measure_Type::CELL_LENGTH_MEASURE,
                            mtk::Primary_Void::PRIMARY,
                            mtk::Master_Slave::MASTER );

            // property type for the SP
            enum class Property_Type
            {
                MATERIAL,
                SLIPLENGTH,    // slip length
                MAX_ENUM
            };

            /*
             *      mParameters( 0 )(0) - alpha_tang
             */

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Dirichlet_Neumann_Nitsche();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Dirichlet_Neumann_Nitsche(){};

            //------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aDofStrings list of strings describing the dof types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dof_type_list(
                    moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                    moris::Cell< std::string >&                  aDofStrings,
                    mtk::Master_Slave                            aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes   a cell of group of dv types
             * @param[ in ] aDvStrings list of strings describing the dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void
            set_dv_type_list(
                    moris::Cell< moris::Cell< PDV_Type > >& aDvTypes,
                    moris::Cell< std::string >&             aDvStrings,
                    mtk::Master_Slave                       aIsMaster = mtk::Master_Slave::MASTER )
            {
                Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsMaster );
            }

            //------------------------------------------------------------------------------
            /**
             * get cluster measure tuples
             * @param[ in ] aClusterMeasureTuples list of tuples describing the cluster measure types
             */
            moris::Cell< std::tuple<
                    fem::Measure_Type,
                    mtk::Primary_Void,
                    mtk::Master_Slave > >
            get_cluster_measure_tuple_list();

            //------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter values
             *
             * mPPVal(0) : gamma_t1 -  stability coeff
             * mPPVal(1) : gamma_t2 -  adjoint coeff
             */
            void eval_SP();

            //------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type >& aDofTypes );

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            void
            eval_dSPdMasterDV( const moris::Cell< PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, "SP_Dirichlet_Neumann_Nitsche::eval_dSPdMasterDV - not implemented." );
            }

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_Dirichlet_Neumann_Nitsche_HPP_ */
