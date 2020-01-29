/*
 * cl_FEM_SP_Ghost_Displacement.hpp
 *
 *  Created on: Nov 15, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_
#define SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_

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

        class SP_Ghost_Displacement : public Stabilization_Parameter
        {

//------------------------------------------------------------------------------
        public:

            // cluster measures
            // FIXME add enum for child class to select the needed ones
            real mMasterVolume     = 0.5; // volume on master
            real mSlaveVolume      = 0.5; // volume on slave
            real mInterfaceSurface = 1.0; // surface on master/slave interface
            real mElementSize      = 1.0; // element size


            enum class SP_Property_Type
            {
                MATERIAL,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, SP_Property_Type > mPropertyMap;

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
             *      mParameters( 1 ) - interpolation order
             */
            SP_Ghost_Displacement();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Ghost_Displacement(){};
//------------------------------------------------------------------------------

            void reset_cluster_measures()
            {
                // FIXME cluster measures to reset volume, surface, ...
                mMasterVolume     = mCluster->compute_cluster_cell_measure(mtk::Primary_Void::INTERP,mtk::Master_Slave::MASTER);
//                mSlaveVolume      = mCluster->compute_cluster_cell_measure(mtk::Primary_Void::INTERP,mtk::Master_Slave::SLAVE);
//                mInterfaceSurface = mCluster->compute_cluster_cell_side_measure(mtk::Primary_Void::PRIMARY,mtk::Master_Slave::SLAVE);
                mElementSize      = std::pow(mMasterVolume/M_PI,0.3333333333) ; // mCluster->compute_element_size( arg on how to compute ) //FIXME: decide what to do here and add the formula to cell info length measure.

                std::cout<<"Element Size = "<<mElementSize<<std::endl;
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
                // FIXME check that property type makes sense?

                // set the property in the property cell
                this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            }

//------------------------------------------------------------------------------
            /**
             * set constitutive model
             * @param[ in ] aConstitutiveModel  a constitutive model pointer
             * @param[ in ] aConstitutiveString a string defining the constitutive model
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                                         std::string                           aConstitutiveString,
                                         mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
            {
                // FIXME check that constitutive string makes sense?

                // set the constitutive model in the constitutive model cell
                this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
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
            void eval_dSPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "SP_Ghost_Displacement::eval_dSPdMasterDV: not implemented." );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_ */
