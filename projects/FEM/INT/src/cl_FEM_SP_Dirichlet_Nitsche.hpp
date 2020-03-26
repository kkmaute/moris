/*
 * cl_FEM_SP_Dirichlet_Nitsche.hpp
 *
 *  Created on: Oct 21, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_

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

        class SP_Dirichlet_Nitsche : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        private:

            // cluster measures for the SP
            real mElementSize = 1.0;

        public:

            // Property type for the SP
            enum class SP_Property_Type
            {
                MATERIAL,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, SP_Property_Type > mPropertyMap;

            // CM type for the SP
            enum class SP_Constitutive_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, SP_Constitutive_Type > mConstitutiveMap;

//------------------------------------------------------------------------------
            /*
             * constructor
             * Rem: mParameters( 0 ) - gamma penalty parameter
             */
            SP_Dirichlet_Nitsche();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Dirichlet_Nitsche(){};

//------------------------------------------------------------------------------
            /**
             * reset the cluster measures required for this SP
             */
            void reset_cluster_measures()
            {
                // evaluate element size from the cluster
                mElementSize = mCluster->compute_cluster_cell_length_measure( mtk::Primary_Void::PRIMARY,
                                                                              mtk::Master_Slave::MASTER );
            }

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
                             "SP_Dirichlet_Nitsche::set_property - Unknown aPropertyString." );

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
             * dPPdMasterDOF ( 1 x numDerDof )
             */
            void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a master dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dPPdMasterDV ( 1 x numDerDv )
             */
            void eval_dSPdMasterDV( const moris::Cell< GEN_DV > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Dirichlet_Nitsche - eval_dSPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_ */
