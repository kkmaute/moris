/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Ghost.hpp
 *
 *  Created on: Jan 21, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_GHOST_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_GHOST_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Ghost : public IWG
        {

            // order of Shape functions
            uint mOrder;

//------------------------------------------------------------------------------
        public:

            enum class IWG_Property_Type
            {
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

            enum class IWG_Constitutive_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

            enum class IWG_Stabilization_Type
            {
                GHOST_DISPL_1,
                GHOST_DISPL_2,
                GHOST_DISPL_3,
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Isotropic_Struc_Linear_Ghost();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Struc_Linear_Ghost(){};

//------------------------------------------------------------------------------
            /**
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                              std::string                                aStabilizationString )
            {
                // check that aConstitutiveString makes sense
                MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                             "IWG_Isotropic_Struc_Linear_Ghost::set_stabilization_parameter - Unknown aStabilizationString." );

                // set the stabilization parameter in the stabilization parameter cell
                this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
            }

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aResidual cell of residual vectors to fill
             */
            void compute_residual( real tWStar );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             */
            void compute_jacobian( real tWStar );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             * @param[ in ] aResidual  cell of residual vectors to fill
             */
            void compute_jacobian_and_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );

//------------------------------------------------------------------------------
            /**
             * method to assemble "normal matrix" from normal vector needed for
             * 2nd and 3rd order Ghost formulations
             * @param[ in ] aOrderGhost order of derivatives and ghost formulation
             * @param[ in ] aNormal     normal matrix to fill
             */
            void get_normal_matrix ( uint aOrderGhost,
                                     Matrix< DDRMat > & aNormal );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */



#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_GHOST_HPP_ */
