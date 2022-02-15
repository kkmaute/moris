/*
 * cl_FEM_SP_Spalart_Allmaras_Turbulence.hpp
 *
 *  Created on: May 19, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_SUPG_Spalart_Allmaras_Turbulence : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private :

                // default dof type
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;
                MSI::Dof_Type mMasterDofVelocity  = MSI::Dof_Type::VX;

                // Spalart-Allmaras model constants
                real mCb2 = 0.6220;
                real mSigma = 2.0/3.0;

                // internal threshold
                const real mEpsilon = MORIS_REAL_EPS;

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        SPALART_ALLMARAS_TURBULENCE,
                        MAX_ENUM
                };

                // flag for evaluation
                bool mLengthScaleEval = true;
                moris::Matrix< DDBMat > mdLengthScaledMasterDofEval;

                // storage
                real mLengthScale;
                moris::Cell< Matrix< DDRMat > > mdLengthScaledMasterDof;

                /*
                 * Rem: mParameters( 0 ) - No parameter used
                 */

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_SUPG_Spalart_Allmaras_Turbulence();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_SUPG_Spalart_Allmaras_Turbulence(){};

                //------------------------------------------------------------------------------
                /**
                 * set function pointers for evaluation
                 */
                void set_function_pointers();

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 * child implementation
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * create a global dof type list including constitutive and property dependencies
                 * child implementation
                 */
                void build_global_dof_type_list();

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                        moris::Cell< std::string >             & aDvStrings,
                        mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsMaster );
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
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdMasterDV - not implemented." );
                }

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------
                /**
                 * return the length scale
                 * @param[ out ] mLengthScale length scale for SUPG
                 */
                real length_scale();

                /**
                 * evaluate the length scale parameter
                 */
                void eval_length_scale();

                //------------------------------------------------------------------------------
                /**
                 * return the length scale derivative wrt to a master dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                const Matrix< DDRMat > & dlengthscaledmasteru( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * evaluate the length scale derivative wrt to a master dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_dlengthscaledmasteru( const moris::Cell< MSI::Dof_Type > & aDofTypes );

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_SPALART_ALLMARAS_TURBULENCE_HPP_ */
