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
             * evaluate the constitutive model stress
             * @param[ in ] aStress a matrix to fill with evaluation
             */
            void eval_stress( Matrix< DDRMat > & aStress );

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
             * evaluate the constitutive model stress derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adStressdDOF a matrix to fill with derivative evaluation
             */
            void eval_dStressdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                   Matrix< DDRMat >             & adStressdDOF );

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
