/*
 * cl_FEM_Constitutive_Model.cpp
 *
 *  Created on: Dec 10, 2019
 *      Author: noel
 */

#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void Constitutive_Model::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager )
        {
            // set the field interpolator manager for the constitutive model
            mFIManager = aFieldInterpolatorManager;

//            // FIXME
//            // get the list of dof types for the CM
//            moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes
//            = this->get_global_dof_type_list();
//
//            // get the number of dof type for the CM
//            uint tNumDofTypes = tCMDofTypes.size();
//
//            // set the size of the field interpolators list for the CM
//            mDofFI.resize( tNumDofTypes, nullptr );
//
//            // loop over the dof types
//            for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//            {
//                // fill the field interpolators list for the CM
//                mDofFI( iDof ) = mFIManager->get_field_interpolators_for_type( tCMDofTypes( iDof )( 0 ) );
//            }
//            // END FIXME
//
//            // FIXME
//            // get the list of dv types for the CM
//            moris::Cell< moris::Cell< MSI::Dv_Type > > tCMDvTypes
//            = this->get_global_dv_type_list();
//
//            // get the number of dv type for the CM
//            uint tNumDvTypes = tCMDvTypes.size();
//
//            // set the size of the field interpolators list for the CM
//            mDvFI.resize( tNumDvTypes, nullptr );
//
//            // loop over the dof types
//            for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
//            {
//                // fill the field interpolator list for the CM
//                mDvFI( iDv ) = mFIManager->get_field_interpolators_for_type( tCMDvTypes( iDv )( 0 ) );
//            }
//            // END FIXME

            // loop over the underlying properties
            for( std::shared_ptr< Property > tProp : this->get_properties() )
            {
                if (tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( mFIManager );
                }
            }
        }

//------------------------------------------------------------------------------
        void Constitutive_Model::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                  moris::Cell< MSI::Dv_Type >  & aDvTypes )
        {
            // init dof counter
            uint tDofCounter = 0;
            uint tDvCounter  = 0;

            // loop over direct dof dependencies
            for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // update counter
                tDofCounter += mDofTypes( iDof ).size();
            }

            // loop over direct dv dependencies
            for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // update counter
                tDvCounter += mDvTypes( iDv ).size();
            }

            // loop over properties
            for ( std::shared_ptr< Property > tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // update counter
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // reserve memory for the non unique dof and dv types
            aDofTypes.reserve( tDofCounter );
            aDvTypes.reserve( tDvCounter );

            // loop over direct dof dependencies
            for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // populate the dof type list
                aDofTypes.append( mDofTypes( iDof ) );
            }

            // loop over direct dv dependencies
            for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // populate the dv type list
                aDvTypes.append( mDvTypes( iDv ) );
            }

            // loop over the properties
            for ( std::shared_ptr< Property > tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof and dv type lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }
        }

//------------------------------------------------------------------------------
            void Constitutive_Model::eval_dFluxdDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                              Matrix< DDRMat >             & adFluxdDOF_FD,
                                                              real                           aPerturbation )
            {
                // get the field interpolator for type
                Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adFluxdDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init dof counter
                uint tDofCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Plus = this->flux();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Minus = this->flux();

                        // evaluate Jacobian
                        adFluxdDOF_FD.get_column( tDofCounter ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
        void Constitutive_Model::eval_dStraindDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                            Matrix< DDRMat >             & adStraindDOF_FD,
                                                            real                           aPerturbation )
        {
            // get the field interpolator for type
            Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of master dofs wrt which derivative is computed
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // FIXME works only for diffusion
            // set size for derivative
            adStraindDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // init dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                    // setting the perturbed coefficients
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mProperties.size();
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // reset constitutive model
                    this->reset_eval_flags();

                    // evaluate the residual
                    Matrix< DDRMat > tStrain_Plus = this->strain();

                    // perturbation of the coefficent
                    tCoeffPert = tCoeff;
                    tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                    // setting the perturbed coefficients
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // reset constitutive model
                    this->reset_eval_flags();

                    // evaluate the residual
                    Matrix< DDRMat > tStrain_Minus = this->strain();

                    // evaluate Jacobian
                    adStraindDOF_FD.get_column( tDofCounter ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

//------------------------------------------------------------------------------
            void Constitutive_Model::eval_dFluxdDV_FD( const moris::Cell< MSI::Dv_Type > & aDvTypes,
                                                             Matrix< DDRMat >            & adFluxdDV_FD,
                                                             real                          aPerturbation )
            {
                // get the field interpolator for type
                Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDvTypes( 0 ) );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDv     = tFI->get_number_of_space_time_coefficients();
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adFluxdDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init dv counter
                uint tDvCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Plus = this->flux();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Minus = this->flux();

                        // evaluate Jacobian
                        adFluxdDV_FD.get_column( tDvCounter ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dv counter
                        tDvCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
            void Constitutive_Model::eval_dStraindDV_FD( const moris::Cell< MSI::Dv_Type > & aDvTypes,
                                                               Matrix< DDRMat >            & adStraindDV_FD,
                                                               real                          aPerturbation )
            {
                // get the field interpolator for type
                Field_Interpolator* tFI = mFIManager->get_field_interpolators_for_type( aDvTypes( 0 ) );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDv     = tFI->get_number_of_space_time_coefficients();
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adStraindDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init dv counter
                uint tDvCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Plus = this->strain();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Minus = this->strain();

                        // evaluate Jacobian
                        adStraindDV_FD.get_column( tDvCounter ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dv counter
                        tDvCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

    }/* end_fem_namespace */
}/* end_moris_namespace */
