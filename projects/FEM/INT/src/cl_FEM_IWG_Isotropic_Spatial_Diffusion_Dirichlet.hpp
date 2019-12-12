/*
 * cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp
 *
 *  Created on: Mar 22, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Isotropic_Spatial_Diffusion_Dirichlet : public IWG
        {

//------------------------------------------------------------------------------
        public:

            enum class IWG_Property_Type
            {
                DIRICHLET,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

            enum class IWG_Constitutive_Type
            {
                DIFF_LIN_ISO,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

            enum class IWG_Stabilization_Type
            {
                DIRICHLET_NITSCHE,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

//------------------------------------------------------------------------------
            /*
             * constructor
             */
            IWG_Isotropic_Spatial_Diffusion_Dirichlet()
            {
                // set size for the property pointer cell
                mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

                // populate the property map
                mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;

                // set size for the constitutive model pointer cell
                mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

                // populate the constitutive map
                mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;

                // set size for the stabilization parameter pointer cell
                mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

                // populate the stabilization map
                mStabilizationMap[ "DirichletNitsche" ] = IWG_Stabilization_Type::DIRICHLET_NITSCHE;
            };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Spatial_Diffusion_Dirichlet(){};

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
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                              std::string                                aStabilizationString )
            {
                // FIXME check that stabilization string makes sense?

                // set the stabilization parameter in the stabilization parameter cell
                this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
            }

//------------------------------------------------------------------------------
            /**
             * computes the residual
             * @param[ in ] aResidual cell of residual vectors to fill
             */
            void compute_residual( real tWStar );

//------------------------------------------------------------------------------
            /**
             * computes the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             */
            void compute_jacobian( real tWStar );

//------------------------------------------------------------------------------
            /**
             * computes the residual and the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             * @param[ in ] aResidual  cell of residual vectors to fill
             */
            void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >                & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_ */
