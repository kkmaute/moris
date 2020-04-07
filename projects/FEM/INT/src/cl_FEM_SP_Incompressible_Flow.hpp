/*
 * cl_FEM_SP_Incompressible_Flow.hpp
 *
 *  Created on: Mar 03, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_
#define SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class SP_Incompressible_Flow : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        public:

            // property type for the SP
            enum class Property_Type
            {
                DENSITY, // fluid density
                VISCOSITY,  // fluid viscosity
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, Property_Type > mPropertyMap;

            // space dimension FIXME
            uint mSpaceDim = 2;

            // pointer to function for G evaluation
            void ( * mEvalGFunc )( Matrix< DDRMat >   & aG,
                                   Matrix< DDRMat >   & aInvSpaceJacobian );

            /*
            * Rem: mParameters( 0 ) - CI = 36
            */

//------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Incompressible_Flow()
            {
                // set the property pointer cell size
                mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

                // populate the map
                mPropertyMap[ "Density" ]   = Property_Type::DENSITY;
                mPropertyMap[ "Viscosity" ] = Property_Type::VISCOSITY;

                // populate the dof map (default)
                mMasterDofMap[ "Velocity" ] = MSI::Dof_Type::VX;
                mMasterDofMap[ "Pressure" ] = MSI::Dof_Type::P;
            }

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Incompressible_Flow(){};

//------------------------------------------------------------------------------
            /**
             * reset the cluster measures required for this SP
             */
            void reset_cluster_measures()
            {
//                // evaluate element size from the cluster
//                mElementSize = mCluster->compute_cluster_cell_length_measure( mtk::Primary_Void::PRIMARY,
//                                                                              mtk::Master_Slave::MASTER );
            }

//------------------------------------------------------------------------------
            /**
             * set function pointers for evaluation
             */
            void set_function_pointers();

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty       a property pointer
             * @param[ in ] aPropertyString a string defining the property
             * @param[ in ] aIsMaster       an enum for master or slave
             */
            void set_property( std::shared_ptr< Property > aProperty,
                               std::string                 aPropertyString,
                               mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check that aPropertyString makes sense
                MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                             "SP_Incompressible_Flow::set_property - Unknown aPropertyString." );

                // set the property in the property cell
                this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            void eval_SP();

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            void eval_dSPdMasterDV( const moris::Cell< GEN_DV > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Incompressible_Flow::eval_dSPdMasterDV - not implemented." );
            }

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
            /**
             * evaluate G
             * where Gij = sum_d dxi_d/dx_i dxi_d/dx_j, d = 1, ..., nSpaceDim
             * @param[ in ] aG a matrix to fill with G
             */
            void eval_G( Matrix< DDRMat > & aG );

            /**
             * evaluate G in 2d
             * @param[ in ] aG                a matrix to fill with G
             * @param[ in ] aInvSpaceJacobian inverse jacobian matrix
             */
            static void eval_G_2d( Matrix< DDRMat > & aG,
                                   Matrix< DDRMat > & aInvSpaceJacobian );

            /**
             * evaluate G in 3d
             * @param[ in ] aG                a matrix to fill with G
             * @param[ in ] aInvSpaceJacobian inverse jacobian matrix
             */
            static void eval_G_3d( Matrix< DDRMat > & aG,
                                   Matrix< DDRMat > & aSpaceJacobian );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_ */
