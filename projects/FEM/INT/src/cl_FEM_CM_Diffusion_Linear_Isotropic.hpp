/*
 * cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class CM_Diffusion_Linear_Isotropic : public Constitutive_Model
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             * trivial constructor
             */
            CM_Diffusion_Linear_Isotropic()
            {
                // set the constitutive type
                mConstitutiveType = fem::Constitutive_Type::DIFF_LIN_ISO;
            };

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Diffusion_Linear_Isotropic(){};

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux
             * @param[ in ] aFlux a matrix to fill with evaluation
             */
            void eval_flux( Matrix< DDRMat > & aFlux );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain
             * @param[ in ] aStrain a matrix to fill with evaluation
             */
            void eval_strain( Matrix< DDRMat > & aStrain );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model matrix
             * @param[ in ] aConst a matrix to fill with evaluation
             */
            void eval_const( Matrix< DDRMat > & aConst );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF a matrix to fill with derivative evaluation
             */
            void eval_dFluxdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                 Matrix< DDRMat >             & adFluxdDOF );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
             */
            void eval_dStraindDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                   Matrix< DDRMat >             & adStraindDOF );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes   a dof type wrt which the derivative is evaluated
             * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
             */
            void eval_dConstdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                   Matrix< DDRMat >            & adConstdDOF );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_ */
