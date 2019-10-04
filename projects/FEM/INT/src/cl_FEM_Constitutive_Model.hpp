/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Constitutive model
         */
        class Constitutive_Model
        {

        protected :

            // constitutive model type
            fem::Constitutive_Type mConstitutiveType;

            // dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // dof type map
            Matrix< DDSMat > mDofTypeMap;

            // field interpolators
            moris::Cell< Field_Interpolator* > mFieldInterpolators;

            // property type list
            moris::Cell< fem::Property_Type > mPropTypes;

            // properties
            moris::Cell< Property* > mProperties;

            // global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mGlobalDofTypes;

            // global dof type map
            Matrix< DDSMat > mGlobalDofTypeMap;

            // spatial dimensions
            uint mSpaceDim;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Constitutive_Model(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Constitutive_Model(){};

//------------------------------------------------------------------------------
            /**
             * returns a constitutive type for the constitutive model
             */
            const fem::Constitutive_Type & get_constitutive_type() const
            {
                // return constitutive type
                return mConstitutiveType;
            };

//------------------------------------------------------------------------------
            /**
             * set space dimension
             * @param[ in ] aSpaceDim a spatial dimension
             */
            void set_space_dim( uint aSpaceDim )
            {
                MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4, "Constitutive_Model::set_space_dim - wrong space dimension.");

                mSpaceDim = aSpaceDim;
            }

//------------------------------------------------------------------------------
            /**
             * sets constitutive model active dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             */
            void set_dof_type_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
            {
                mDofTypes = aDofTypes;

                this->build_dof_type_map();
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types active for the constitutive model
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list() const
            {
                // return dof type list
                return mDofTypes;
            };

//------------------------------------------------------------------------------
            void build_dof_type_map()
            {
                // determine the max Dof_Type enum
                sint tMaxEnum = 0;
                for( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDof )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dof_Type map size
                mDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dof_Type map
                for( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
                {
                    // fill the property map
                    mDofTypeMap( static_cast< int >( mDofTypes( iDof )( 0 ) ), 0 ) = iDof;
                }
            }

//------------------------------------------------------------------------------
            /**
             * sets constitutive model active property types
             * @param[ in ] aPropertyTypes a cell of cell of property types
             *
             */
            void set_property_type_list( const moris::Cell< fem::Property_Type  > & aPropertyTypes )
            {
                mPropTypes = aPropertyTypes;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of property type active for the constitutive model
             */
            const moris::Cell< fem::Property_Type > & get_property_type_list() const
            {
                // return property type list
                return mPropTypes;
            };

//------------------------------------------------------------------------------
            /**
             * sets field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             */
            void set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == mGlobalDofTypes.size(),
                              "Constitutive_Model::set_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == mGlobalDofTypes( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Constitutive_Model::set_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                mFieldInterpolators = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * gets field interpolators
             */
            const moris::Cell< Field_Interpolator* > & get_field_interpolators()
            {
                // return field interpolator pointers
                return mFieldInterpolators;
            };

//------------------------------------------------------------------------------
            /**
             * checks that field interpolators were assigned
             */
            void check_field_interpolators()
            {
                // check field interpolators cell size
                MORIS_ASSERT( mFieldInterpolators.size() == mGlobalDofTypes.size(),
                              "Constitutive_Model::check_field_interpolators - wrong FI size. " );

               // loop over the field interpolator pointers
               for( uint iFI = 0; iFI < mGlobalDofTypes.size(); iFI++ )
               {
                   // check that the field interpolator was set
                   MORIS_ASSERT( mFieldInterpolators( iFI ) != nullptr,
                                 "Constitutive_Model::check_field_interpolators - FI missing. " );
               }
            }

//------------------------------------------------------------------------------
            /**
             * sets properties
             * @param[ in ] aProperties cell of property pointers
             */
            void set_properties( moris::Cell< Property* > & aProperties )
            {
                // check input size
                MORIS_ASSERT( aProperties.size() == mPropTypes.size(),
                              "Constitutive_Model::set_properties - master, wrong input size. " );

                // check property type
                bool tCheckProp = true;
                for( uint iProp = 0; iProp < aProperties.size(); iProp++ )
                {
                    tCheckProp = tCheckProp && ( aProperties( iProp )->get_property_type() == mPropTypes( iProp ) );
                }
                MORIS_ASSERT( tCheckProp, "Constitutive_Model::set_properties - wrong property type. ");

                // set properties
                mProperties = aProperties;

                // create a global dof type list
                this->build_global_dof_type_list();
            }

//------------------------------------------------------------------------------
            /**
             * gets properties
             */
            const moris::Cell< Property* > & get_properties()
            {
                // return property pointers
                return mProperties;
            };

//------------------------------------------------------------------------------
            /**
             * checks that properties were assigned
             */
            void check_properties()
            {
                // check master properties cell size
                MORIS_ASSERT( mProperties.size() == mPropTypes.size(),
                              "Constitutive_Model::check_properties - wrong property size. " );

                // loop over all master properties and check that they are assigned
                for( uint iProp = 0; iProp < mPropTypes.size(); iProp++ )
                {
                    MORIS_ASSERT( mProperties( iProp ) != nullptr,
                                  "Constitutive_Model::check_properties - property missing. " );
                }
            }

//------------------------------------------------------------------------------
            /**
             * creates a global dof type list including constitutive model and properties dependencies
             */
            void build_global_dof_type_list()
            {
                // set the size of the dof type list
                uint tCounterMax = mDofTypes.size();

                for ( Property* tProperty : mProperties )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }
                mGlobalDofTypes.resize( tCounterMax );
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dof counter
                uint tCounter = 0;

                // get active dof type for IWG
                for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
                {
                    tCheckList( tCounter ) = static_cast< uint >( mDofTypes( iDOF )( 0 ) );
                    mGlobalDofTypes( tCounter ) = mDofTypes( iDOF );
                    tCounter++;
                }

                for ( Property* tProperty : mProperties )
                {
                    // get active dof type
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
                            mGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );
                            tCounter++;
                        }
                    }
                }

                // get the number of unique dof type groups, i.e. the number of interpolators
                mGlobalDofTypes.resize( tCounter );

                // build global dof type map
                this->build_global_dof_type_map();

            };

//------------------------------------------------------------------------------
            /**
             * gets global dof type list
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list()
            {
                // return global dof type list
                return mGlobalDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * builds global dof type list map
             */
            void build_global_dof_type_map()
            {
                // get number of global dof types
                uint tNumDofTypes = mGlobalDofTypes.size();

                // determine the max Dof_Type enum
                sint tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dof_Type map size
                mGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dof_Type map
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    // fill the property map
                    mGlobalDofTypeMap( static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
                }
            }


//------------------------------------------------------------------------------
            /**
             * checks dependency on a given group of dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             *
             */
            bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // set bool for dependency
                bool tDofDependency = false;

                // if aDofType is an active dof type for the property
                if( static_cast< uint >( aDofType( 0 ) ) < mGlobalDofTypeMap.numel() && mGlobalDofTypeMap(static_cast< uint >( aDofType( 0 ) ) ) != -1 )
                {
                    // bool is set to true
                    tDofDependency = true;
                }
                // return bool for dependency
                return tDofDependency;
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model flux
             * @param[ in ] aFlux a matrix to fill with flux evaluation
             */
            virtual void eval_flux( Matrix< DDRMat > & aFlux )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_flux - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model strain
             * @param[ in ] aStrain a matrix to fill with strain evaluation
             */
            virtual void eval_strain( Matrix< DDRMat > & aStrain )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_strain - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model constitutive matrix
             * @param[ in ] aConst a matrix to fill with constitutive matrix evaluation
             */
            virtual void eval_const( Matrix< DDRMat > & aConst )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_const - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dFluxdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                         Matrix< DDRMat >             & adFluxdDOF )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dStraindDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                           Matrix< DDRMat >             & adStraindDOF )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model constitutive matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dConstdDOF( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                          Matrix< DDRMat >             & adConstdDOF )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model stress derivative wrt to a dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dFluxdDOF_FD( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                      Matrix< DDRMat >           & adFluxdDOF_FD,
                                      real                         aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of master dofs wrt which derivative is computed
                uint tDerNumDof = mFieldInterpolators( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adFluxdDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mFieldInterpolators( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDof; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mFieldInterpolators( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // evaluate the residual
                    Matrix< DDRMat > tFlux_Plus;
                    this->eval_flux( tFlux_Plus );

                    // perturbation of the coefficent
                    tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mFieldInterpolators( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // evaluate the residual
                    Matrix< DDRMat > tFlux_Minus;
                    this->eval_flux( tFlux_Minus );

                    // evaluate Jacobian
                    adFluxdDOF_FD.get_column( iCoeff ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
                }
                // reset the coefficients values
                mFieldInterpolators( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
            /**
            * evaluates the constitutive model strain derivative wrt to a dof type
            * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
            * @param[ in ] adStraindDOF_FD a matrix to fill with derivative evaluation
            * @param[ in ] aPerturbation   real to perturb for FD
            */
            void eval_dStraindDOF_FD( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                      Matrix< DDRMat >             & adStraindDOF_FD,
                                      real                           aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of master dofs wrt which derivative is computed
                uint tDerNumDof = mFieldInterpolators( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adStraindDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mFieldInterpolators( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDof; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mFieldInterpolators( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // evaluate the residual
                    Matrix< DDRMat > tStrain_Plus;
                    this->eval_strain( tStrain_Plus );

                    // perturbation of the coefficent
                    tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mFieldInterpolators( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                    {
                        mProperties( iProp )->reset_eval_flags();
                    }

                    // evaluate the residual
                    Matrix< DDRMat > tStrain_Minus;
                    this->eval_strain( tStrain_Minus );

                    // evaluate Jacobian
                    adStraindDOF_FD.get_column( iCoeff ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
                }
                // reset the coefficients values
                mFieldInterpolators( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_ */
