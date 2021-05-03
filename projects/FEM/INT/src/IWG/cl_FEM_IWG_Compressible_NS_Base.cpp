/*
 * cl_FEM_IWG_Compressible_NS_Base.cpp
 *
 *  Created on: Apr 23, 2021
 *      Author: wunsch
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Base.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void IWG_Compressible_NS_Base::reset_spec_eval_flags()
        {
            // reset eval flags
            mYEval = true;
            mdYdtEval = true;
            mdYdxEval = true;
            md2Ydx2Eval = true;

            mWEval = true;
            mdWdtEval = true;
            mdWdxEval = true;
            md2Wdx2Eval = true;

            mAEval = true;
            mKEval = true;
            mKijiEval = true;

            // reset flags for child
            this->reset_child_eval_flags();
        }

        //------------------------------------------------------------------------------

        uint IWG_Compressible_NS_Base::num_space_dims()
        {
            // get number of spatial dimensions from velocity, momentum, etc. field interpolator
            return mMasterFIManager->get_field_interpolators_for_type( mVectorDof )->get_number_of_fields();
        }

        //------------------------------------------------------------------------------

        uint IWG_Compressible_NS_Base::num_bases()
        {
            // get number of bases from first state var FI
            return mMasterFIManager->get_field_interpolators_for_type( mFirstDof )->get_number_of_space_time_bases();
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        
        const Matrix< DDRMat > & IWG_Compressible_NS_Base::Y()
        {
            // check if the variable vectors have already been assembled
            if( !mYEval )
            {      
                return mY;
            }   
            mYEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Base::Y() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator * tFI1 =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof  );
            Field_Interpolator * tFI2 =  mMasterFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator * tFI3 =  mMasterFIManager->get_field_interpolators_for_type( mLastDof   );

            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    mY = {  { tFI1->val()( 0 ) }, 
                            { tFI2->val()( 0 ) },
                            { tFI2->val()( 1 ) }, 
                            { tFI3->val()( 0 ) } };
                    break;
                }

                // 3D
                case 3 :
                {
                    mY = {  { tFI1->val()( 0 ) }, 
                            { tFI2->val()( 0 ) },
                            { tFI2->val()( 1 ) },
                            { tFI2->val()( 2 ) }, 
                            { tFI3->val()( 0 ) } };
                    break;
                }
            
                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::Y() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }

            // return Y
            return mY;
        }

        //------------------------------------------------------------------------------
        
        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dYdt()
        {
            // check if the variable vectors have already been assembled
            if( !mdYdtEval )
            {      
                return mdYdt;
            }   
            mdYdtEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Base::dYdt() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator * tFI1 =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof  );
            Field_Interpolator * tFI2 =  mMasterFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator * tFI3 =  mMasterFIManager->get_field_interpolators_for_type( mLastDof   );

            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    mdYdt = {  
                            { tFI1->gradt( 1 )( 0 ) }, 
                            { tFI2->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 1 ) }, 
                            { tFI3->gradt( 1 )( 0 ) } };
                    break;
                }

                // 3D
                case 3 :
                {
                    mdYdt = {  
                            { tFI1->gradt( 1 )( 0 ) }, 
                            { tFI2->gradt( 1 )( 0 ) },
                            { tFI2->gradt( 1 )( 1 ) },
                            { tFI2->gradt( 1 )( 2 ) }, 
                            { tFI3->gradt( 1 )( 0 ) } };
                    break;
                }
            
                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::dYdt() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }

            // return dYdt
            return mdYdt;
        }

        //------------------------------------------------------------------------------
        
        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dYdx( const uint aSpatialDirection )
        {
            // check if the variable vector has already been assembled
            if( !mdYdxEval )
            {      
                return mdYdx( aSpatialDirection );
            }  
            mdYdxEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Base::dYdx() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator * tFI1 =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof  );
            Field_Interpolator * tFI2 =  mMasterFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator * tFI3 =  mMasterFIManager->get_field_interpolators_for_type( mLastDof   );

            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    // set size of cell
                    mdYdx.resize( 2 );

                    // x-derivative
                    mdYdx( 0 ) = {  
                            { tFI1->gradx( 1 )( 0 )    }, 
                            { tFI2->gradx( 1 )( 0, 0 ) },
                            { tFI2->gradx( 1 )( 0, 1 ) }, 
                            { tFI3->gradx( 1 )( 0 )    } };
 
                    // y-derivative
                    mdYdx( 1 ) = {  
                            { tFI1->gradx( 1 )( 1 )    }, 
                            { tFI2->gradx( 1 )( 1, 0 ) },
                            { tFI2->gradx( 1 )( 1, 1 ) }, 
                            { tFI3->gradx( 1 )( 1 )    } };

                    break;
                }

                // 3D
                case 3 :
                {
                    // set size of cell
                    mdYdx.resize( 3 );

                    // x-derivative
                    mdYdx( 0 ) = {  
                            { tFI1->gradx( 1 )( 0 )    }, 
                            { tFI2->gradx( 1 )( 0, 0 ) },
                            { tFI2->gradx( 1 )( 0, 1 ) }, 
                            { tFI2->gradx( 1 )( 0, 2 ) }, 
                            { tFI3->gradx( 1 )( 0 )    } };
 
                    // y-derivative
                    mdYdx( 1 ) = {  
                            { tFI1->gradx( 1 )( 1 )    }, 
                            { tFI2->gradx( 1 )( 1, 0 ) },
                            { tFI2->gradx( 1 )( 1, 1 ) }, 
                            { tFI2->gradx( 1 )( 1, 2 ) }, 
                            { tFI3->gradx( 1 )( 1 )    } };

                    // z-derivative
                    mdYdx( 2 ) = {  
                            { tFI1->gradx( 1 )( 2 )    }, 
                            { tFI2->gradx( 1 )( 2, 0 ) },
                            { tFI2->gradx( 1 )( 2, 1 ) }, 
                            { tFI2->gradx( 1 )( 2, 2 ) }, 
                            { tFI3->gradx( 1 )( 2 )    } };     

                    break;
                }
            
                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::dYdx() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }

            // return dYdt
            return mdYdx( aSpatialDirection );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::d2Ydx2( const uint aI, const uint aJ )
        {
            // convert the two indices into one for condensed tensor
            uint tFlatIndex = convert_index_pair_to_flat( aI, aJ, this->num_space_dims() );
            
            // check if the variable vectors have already been assembled
            if( !md2Ydx2Eval )
            {      
                return md2Ydx2( tFlatIndex );
            }   
            md2Ydx2Eval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType ), 
                    "IWG_Compressible_NS_Base::d2Ydx2() - check for residual DoF types failed. See Error message above for more info." );

            // get field interpolators
            Field_Interpolator * tFI1 =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof  );
            Field_Interpolator * tFI2 =  mMasterFIManager->get_field_interpolators_for_type( mVectorDof );
            Field_Interpolator * tFI3 =  mMasterFIManager->get_field_interpolators_for_type( mLastDof   );

            // construct Y - vector based on the number of space dims
            switch ( this->num_space_dims() )
            {
                // 2D
                case 2 :
                {
                    // set size of cell
                    md2Ydx2.resize( 3 );

                    // loop over six derivative combinations
                    for ( uint iDeriv = 0; iDeriv < 3; iDeriv++ )
                    {
                        md2Ydx2( iDeriv ) = {  
                                { tFI1->gradx( 2 )( iDeriv )    }, 
                                { tFI2->gradx( 2 )( iDeriv, 0 ) },
                                { tFI2->gradx( 2 )( iDeriv, 1 ) }, 
                                { tFI3->gradx( 2 )( iDeriv )    } }; 
                    }    
                    break;
                }

                // 3D
                case 3 :
                {
                    // set size of cell
                    md2Ydx2.resize( 6 );

                    // loop over six derivative combinations
                    for ( uint iDeriv = 0; iDeriv < 6; iDeriv++ )
                    {
                        md2Ydx2( iDeriv ) = {  
                                { tFI1->gradx( 2 )( iDeriv )    }, 
                                { tFI2->gradx( 2 )( iDeriv, 0 ) },
                                { tFI2->gradx( 2 )( iDeriv, 1 ) }, 
                                { tFI2->gradx( 2 )( iDeriv, 2 ) }, 
                                { tFI3->gradx( 2 )( iDeriv )    } }; 
                    } 
                    break;
                }
            
                // invalid number of spatial dimensions
                default :
                {
                    MORIS_ERROR( false, "IWG_Compressible_NS_Base::d2Ydx2() - Number of spatial dimensions must be 2 or 3;" );
                    break;
                }
            }

            // return value
            return md2Ydx2( tFlatIndex );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        
        const Matrix< DDRMat > & IWG_Compressible_NS_Base::W()
        {
            // check if the variable vectors have already been assembled
            if( !mWEval )
            {      
                return mW;
            } 
            mWEval = false;  

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::W() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;
            
            // initialize W matrix
            mW.set_size( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );   

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Field_Interpolator * tFI =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof );
            Matrix< DDRMat > tN = tFI->N()( { 0, 0 }, { 0, this->num_bases() - 1 } );
            
            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // put test functions into matrices
                mW( { iVar, iVar  }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) = 
                        tN.matrix_data();
            }

            // return value
            return mW;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdt()
        {
            // check if the variable vectors have already been assembled
            if( !mdWdtEval )
            {      
                return mdWdt;
            }  
            mdWdtEval = false; 

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::dWdt() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;
            
            // initialize W matrix
            mdWdt.set_size( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );   

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Field_Interpolator * tFI =  mMasterFIManager->get_field_interpolators_for_type( mFirstDof );
            Matrix< DDRMat > tdNdt = tFI->dnNdtn( 1 )( { 0, 0 }, { 0, this->num_bases() - 1 } );
            
            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // put test functions into matrices
                mdWdt( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) = 
                        tdNdt.matrix_data(); 
            }

            // return value
            return mdWdt;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::dWdx( const uint aSpatialDirection )
        {
            // check if the variable vectors have already been assembled
            if( !mdWdxEval )
            {      
                return mdWdx( aSpatialDirection );
            }   
            mdWdxEval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::dWdx() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;

            // create zero mat
            Matrix< DDRMat > tZeroMat( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );   
            
            // initialize W matrix
            mdWdx.assign( this->num_space_dims(), tZeroMat );

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Matrix< DDRMat > tdNdx = mMasterFIManager->get_field_interpolators_for_type( mFirstDof )->dnNdxn( 1 );
            
            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
                {
                    mdWdx( iDim )( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) = 
                            tdNdx( { iDim, iDim }, { 0, this->num_bases() - 1 } );
                }
            }

            // return value
            return mdWdx( aSpatialDirection );
        }

        //------------------------------------------------------------------------------
        
        const Matrix< DDRMat > & IWG_Compressible_NS_Base::d2Wdx2( const uint aI, const uint aJ )
        {
            // convert the two indices into one for condensed tensor
            uint tFlatIndex = convert_index_pair_to_flat( aI, aJ, this->num_space_dims() );
            
            // check if the variable vectors have already been assembled
            if( !md2Wdx2Eval )
            {      
                return md2Wdx2( tFlatIndex );
            } 
            md2Wdx2Eval = false;

            // check residual DoF types
            MORIS_ASSERT( check_residual_dof_types( mResidualDofType  ), 
                    "IWG_Compressible_NS_Bulk::d2Wdx2() - check for residual DoF types failed." );

            // get number of state variable fields
            uint tNumStateVars = this->num_space_dims() + 2;
            uint tNumSecondSpaceDerivs = 3 * this->num_space_dims() - 3;
            
            // create zero mat
            Matrix< DDRMat > tZeroMat( tNumStateVars, tNumStateVars * this->num_bases(), 0.0 );   
            
            // initialize W matrix
            md2Wdx2.assign( tNumSecondSpaceDerivs, tZeroMat ); 

            // get representative values for the different basis function vectors
            // NOTE: only works under the assumption that all state variable fields are interpolated on the same mesh
            Matrix< DDRMat > td2Ndx2 = mMasterFIManager->get_field_interpolators_for_type( mFirstDof )->dnNdxn( 2 );
            
            // go through residual dof types and assemble the test function matrix
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                for ( uint iRow = 0; iRow < tNumSecondSpaceDerivs; iRow++ )
                {
                    md2Wdx2( iRow )( { iVar, iVar }, { iVar * this->num_bases(), ( iVar + 1 ) * this->num_bases() - 1 } ) = 
                            td2Ndx2( { iRow, iRow }, { 0, this->num_bases() - 1 } );
                }
            }

            // return value
            return md2Wdx2( tFlatIndex );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::A( const uint aK )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( ( aK >= 0 ) and ( aK <= this->num_space_dims() ), 
                    "IWG_Compressible_NS_Base::A() - index out of bounds." );

            // check if the variable vectors have already been assembled
            if( !mAEval )
            {
                return mA( aK );
            }  
            
            // set the eval flag
            mAEval = false;          

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mMasterMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // evaluate A-matrices and store them
            eval_A( tMM, tCM, mMasterFIManager, mResidualDofType, mA );  

            // return requested value
            return mA( aK );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::K( const uint aI, const uint aJ )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( ( aI >= 0 ) and ( aI < this->num_space_dims() ) and ( aJ >= 0 ) and ( aJ < this->num_space_dims() ), 
                    "IWG_Compressible_NS_Base::K() - indices out of bounds." );

            // check if K matrices have already been evaluated
            if ( !mKEval )
            {
                return mK( aI )( aJ );
            }

            // set the eval flag
            mKEval = false;            

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropThermalConductivity = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // eval K matrices and store them
            eval_K( tPropDynamicViscosity, tPropThermalConductivity, mMasterFIManager, mK );

            // return requested K matrix
            return mK( aI )( aJ );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Base::Kiji( const uint aJ )
        {
            // check that indices are not out of bounds
            MORIS_ASSERT( ( aJ >= 0 ) and ( aJ < this->num_space_dims() ), 
                    "IWG_Compressible_NS_Base::Kiji() - index out of bounds." );

            // check if Kiji matrices have already been evaluated
            if ( !mKijiEval )
            {
                return mKiji( aJ );
            }

            // set the eval flag
            mKijiEval = false;            

            // get the viscosity
            std::shared_ptr< Property > tPropDynamicViscosity = mMasterProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropThermalConductivity = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // eval spatial derivatives of K matrices and store them
            eval_dKijdxi( tPropDynamicViscosity, tPropThermalConductivity, mMasterFIManager, mKiji );

            // return requested Kiji matrix
            return mKiji( aJ );
        } 


        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */