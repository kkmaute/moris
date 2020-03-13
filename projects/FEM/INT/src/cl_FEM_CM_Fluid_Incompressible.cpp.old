
#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_flux()
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_flux - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_traction - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_testTraction( const Matrix< DDRMat > & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_testTraction - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_strain()
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_strain - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_testStrain()
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_testStrain - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_const()
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_const - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dFluxdDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                          const Matrix< DDRMat >             & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dTractiondDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                              const Matrix< DDRMat >             & aNormal,
                                                              const Matrix< DDRMat >             & aJump )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dTestTractiondDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dStraindDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dConstdDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::flatten_normal( const Matrix< DDRMat > & aNormal,
                                                            Matrix< DDRMat > & aFlatNormal )
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    aFlatNormal.set_size( 2, 3, 0.0 );
                    aFlatNormal( 0, 0 ) = aNormal( 0,0 );
                    aFlatNormal( 0, 2 ) = aNormal( 1,0 );
                    aFlatNormal( 1, 1 ) = aNormal( 1,0 );
                    aFlatNormal( 1, 2 ) = aNormal( 0,0 );
                    break;
                }
                case( 3 ):
                {
                    aFlatNormal.set_size( 3, 6, 0.0 );
                    aFlatNormal( 0, 0 ) = aNormal( 0,0 );
                    aFlatNormal( 1, 1 ) = aNormal( 1,0 );
                    aFlatNormal( 2, 2 ) = aNormal( 2,0 );
                    aFlatNormal( 0, 4 ) = aNormal( 2,0 );
                    aFlatNormal( 0, 5 ) = aNormal( 1,0 );
                    aFlatNormal( 1, 3 ) = aNormal( 2,0 );
                    aFlatNormal( 1, 5 ) = aNormal( 0,0 );
                    aFlatNormal( 2, 3 ) = aNormal( 1,0 );
                    aFlatNormal( 2, 4 ) = aNormal( 0,0 );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "CM_Fluid_Incompressible::flatten_normal - Flattening for normal only implemented in 2 and 3D" );
                }
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::set_space_dim( uint aSpaceDim )
        {
            // check that space dimension is 1, 2, 3
            MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4,
                         "Constitutive_Model::set_space_dim - wrong space dimension.");

            // set space dimension
            mSpaceDim = aSpaceDim;
        }

//--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
