/*
 * cl_FEM_Material_Model.hpp
 *
 *  Created on: Feb 1, 2021
 *      Author: wunsch
 */
#ifndef SRC_FEM_CL_FEM_MATERIAL_MODEL_HPP_
#define SRC_FEM_CL_FEM_MATERIAL_MODEL_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//#include "linalg_typedefs.hpp"
//MRS/CON/src
#include "cl_Cell.hpp"
//LNA/src
#include "cl_Matrix.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
//FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
#include <map>
//GEN/src
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace fem
    {
        class Field_Interpolator_Manager;
        class Set;

        //------------------------------------------------------------------------------
        /**
         * -- Thermodynamic Material Model - ("Material_Model") --
         * This class is used to assign a thermodynamic material model to a phase/set.
         * It is able to compute the internal energy based on the first equation of 
         * state for a given material, and the thermodynamic variables (rho,p,T) based 
         * on the fields inputed (usually 2 of 3 primitive variables or entropy vars.)
         * 
         */
        class Material_Model
        {

            protected :

                // field interpolator manager
                Field_Interpolator_Manager * mFIManager = nullptr;

                // fem set pointer
                Set * mSet = nullptr;

                // dof type list
                moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

                // dof type map
                Matrix< DDSMat > mDofTypeMap;

                // global dof type list
                moris::Cell< moris::Cell< MSI::Dof_Type > > mGlobalDofTypes;

                // global dof type map
                Matrix< DDSMat > mGlobalDofTypeMap;

                // properties
                moris::Cell< std::shared_ptr< Property > > mProperties;

                // local string to property enum map
                std::map< std::string, uint > mPropertyMap;

                // spatial dimensions
                uint mSpaceDim;

                // storage for variables (computed using 2nd EOS)
                Matrix< DDRMat > mDensity;
                Matrix< DDRMat > mDensityDot;
                Matrix< DDRMat > mdDensitydx;
                Matrix< DDRMat > md2Densitydx2;

                moris::Cell< Matrix< DDRMat > > mDensityDof;
                moris::Cell< Matrix< DDRMat > > mDensityDotDof;
                moris::Cell< Matrix< DDRMat > > mdDensitydxDof;
                moris::Cell< Matrix< DDRMat > > md2Densitydx2Dof;

                Matrix< DDRMat > mPressure;
                Matrix< DDRMat > mPressureDot;
                Matrix< DDRMat > mdPressuredx;
                Matrix< DDRMat > md2Pressuredx2;

                moris::Cell< Matrix< DDRMat > > mPressureDof;
                moris::Cell< Matrix< DDRMat > > mPressureDotDof;
                moris::Cell< Matrix< DDRMat > > mdPressuredxDof;
                moris::Cell< Matrix< DDRMat > > md2Pressuredx2Dof;

                Matrix< DDRMat > mTemperature;
                Matrix< DDRMat > mTemperatureDot;
                Matrix< DDRMat > mdTemperaturedx;
                Matrix< DDRMat > md2Temperaturedx2;

                moris::Cell< Matrix< DDRMat > > mTemperatureDof;
                moris::Cell< Matrix< DDRMat > > mTemperatureDotDof;
                moris::Cell< Matrix< DDRMat > > mdTemperaturedxDof;
                moris::Cell< Matrix< DDRMat > > md2Temperaturedx2Dof;                                
                             
                // storage for specific internal energy (computed using 1st EOS)
                Matrix< DDRMat > mEint;
                Matrix< DDRMat > mEintDot;
                Matrix< DDRMat > mdEintdx;
                Matrix< DDRMat > md2Eintdx2;

                moris::Cell< Matrix< DDRMat > > mEintDof;
                moris::Cell< Matrix< DDRMat > > mEintDotDof;
                moris::Cell< Matrix< DDRMat > > mdEintdxDof;
                moris::Cell< Matrix< DDRMat > > md2Eintdx2Dof;

                // storage for other material properties
                Matrix< DDRMat > mAlphaP;
                Matrix< DDRMat > mBetaT;
                Matrix< DDRMat > mCv;
                Matrix< DDRMat > mCp;
                Matrix< DDRMat > mGamma;

                moris::Cell< Matrix< DDRMat > > mAlphaPDof;
                moris::Cell< Matrix< DDRMat > > mBetaTDof;

                // FIXME: for now assume these dof derivs are 0
                moris::Cell< Matrix< DDRMat > > mCvDof;
                moris::Cell< Matrix< DDRMat > > mCpDof;
                moris::Cell< Matrix< DDRMat > > mGammaDof;

                // constitutive model name for input and debug
                std::string mName = "Undefined";

            private:

                // bool for global dof type list and map build
                bool mGlobalDofBuild = true;
                bool mGlobalDofMapBuild = true;

                // flags for thermodynamic variables
                bool mDensityEval = true;
                bool mPressureEval = true;
                bool mTemperatureEval = true;

                bool mDensityDotEval = true;
                bool mPressureDotEval = true;
                bool mTemperatureDotEval = true;

                bool mdDensitydxEval = true;
                bool mdPressuredxEval = true;
                bool mdTemperaturedxEval = true;

                bool md2Densitydx2Eval = true;
                bool md2Pressuredx2Eval = true;
                bool md2Temperaturedx2Eval = true;

                moris::Matrix< DDBMat > mDensityDofEval;
                moris::Matrix< DDBMat > mPressureDofEval;
                moris::Matrix< DDBMat > mTemperatureDofEval;

                moris::Matrix< DDBMat > mDensityDotDofEval;
                moris::Matrix< DDBMat > mPressureDotDofEval;
                moris::Matrix< DDBMat > mTemperatureDotDofEval;                

                moris::Matrix< DDBMat > mdDensitydxDofEval;
                moris::Matrix< DDBMat > mdPressuredxDofEval;
                moris::Matrix< DDBMat > mdTemperaturedxDofEval;

                moris::Matrix< DDBMat > md2Densitydx2DofEval;
                moris::Matrix< DDBMat > md2Pressuredx2DofEval;
                moris::Matrix< DDBMat > md2Temperaturedx2DofEval;                

                // flag for specific internal energy (computed using 1st EOS)
                bool mEintEval = true;
                bool mEintDotEval = true;
                bool mdEintdxEval = true;
                bool md2Eintdx2Eval = true;

                moris::Matrix< DDBMat > mEintDofEval;
                moris::Matrix< DDBMat > mEintDotDofEval;
                moris::Matrix< DDBMat > mdEintdxDofEval;
                moris::Matrix< DDBMat > md2Eintdx2DofEval;

                // flags for other variables
                bool mAlphaPEval = true;
                bool mBetaTEval = true;
                bool mCvEval = true;
                bool mCpEval = true;
                bool mGammaEval = true;

                moris::Matrix< DDBMat > mAlphaPDofEval;
                moris::Matrix< DDBMat > mBetaTDofEval;

                // FIXME: for now assume these dof derivs are 0
                moris::Matrix< DDBMat > mCvDofEval;
                moris::Matrix< DDBMat > mCpDofEval;
                moris::Matrix< DDBMat > mGammaDofEval;

                // booleans indicating dependent 
                bool mCvIsDependent = true;
                bool mCpIsDependent = true;
                bool mGammaIsDependent = true;

                // booleans indicating dependent variables
                bool mDensityIsDependent = true;
                bool mPressureIsDependent = true;
                bool mTemperatureIsDependent = true;

                // Dof-types for Thermodynamic variables
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // function pointers for density
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Density )() = 
                    &Material_Model::density_dep;
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_DensityDot )() = 
                    &Material_Model::DensityDot_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnDensitydxn )( uint aOrder ) = 
                    &Material_Model::dnDensitydxn_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_DensityDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::DensityDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_DensityDotDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::DensityDotDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnDensitydxnDof )( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder ) = 
                    &Material_Model::dnDensitydxnDOF_dep;

                // function pointers for pressure
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Pressure )() = 
                    &Material_Model::pressure_dep;
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_PressureDot )() = 
                    &Material_Model::PressureDot_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnPressuredxn )( uint aOrder ) = 
                    &Material_Model::dnPressuredxn_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_PressureDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::PressureDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_PressureDotDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::PressureDotDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnPressuredxnDof )( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder ) = 
                    &Material_Model::dnPressuredxnDOF_dep;

                // function pointers for temperature
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Temperature )() = 
                    &Material_Model::temperature_dep;
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_TemperatureDot )() = 
                    &Material_Model::TemperatureDot_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnTemperaturedxn )( uint aOrder ) = 
                    &Material_Model::dnTemperaturedxn_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_TemperatureDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::TemperatureDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_TemperatureDotDof )( const moris::Cell< MSI::Dof_Type > & aDofType ) = 
                    &Material_Model::TemperatureDotDOF_dep;
                const Matrix< DDRMat > & ( Material_Model:: * m_get_dnTemperaturedxnDof )( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder ) = 
                    &Material_Model::dnTemperaturedxnDOF_dep;  

                // function pointers for cv and cp (FIXME: for now, cv is given, cp and gamma are dependent)
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Cv )() = 
                    &Material_Model::Cv_triv;   
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Cp )() = 
                    &Material_Model::Cp_dep;    
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_Gamma )() = 
                    &Material_Model::Gamma_dep;   

                // function pointers for cv and cp (FIXME: for now, cv is given, cp and gamma are dependent)
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_CvDof )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = 
                    &Material_Model::CvDOF_triv;   
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_CpDof )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = 
                    &Material_Model::CpDOF_dep;    
                const Matrix< DDRMat > &  ( Material_Model:: * m_get_GammaDof )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = 
                    &Material_Model::GammaDOF_dep;                                                                        

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                Material_Model();

                //------------------------------------------------------------------------------
                /**
                 * virtual destructor
                 */
                virtual ~Material_Model(){};

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       property shared pointer to set
                 * @param[ in ] aPropertyString string describing the property to set
                 */
                void set_property(
                        std::shared_ptr< fem::Property > aProperty,
                        std::string                      aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * get property
                 * @param[ in ]  aPropertyString string describing the property to get
                 * @param[ out ] aProperty       property shared pointer
                 */
                std::shared_ptr< fem::Property > & get_property( std::string aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * set local properties
                 */
                virtual void set_local_properties(){}

                //------------------------------------------------------------------------------
                /**
                 * set name
                 * param[ in ] aName a string for MM name
                 */
                void set_name( std::string aName )
                {
                    mName = aName;
                }

                //------------------------------------------------------------------------------
                /**
                 * get name
                 * param[ out ] mName a string for MM name
                 */
                std::string get_name()
                {
                    return mName;
                }

                //------------------------------------------------------------------------------
                /**
                 * print names
                 */
                void print_names();

                //------------------------------------------------------------------------------
                /*
                 * set member set pointer
                 * @param[ in ] aSetPointer a FEM set pointer
                 */
                void set_set_pointer( Set * aSetPointer )
                {
                    mSet = aSetPointer;
                }

                //------------------------------------------------------------------------------
                /**
                 * set space dimension
                 * @param[ in ] aSpaceDim a spatial dimension
                 */
                virtual void set_space_dim( uint aSpaceDim )
                {
                    // check that space dimension is 1, 2, 3
                    MORIS_ERROR(
                            aSpaceDim > 0 && aSpaceDim < 4,
                            "Material_Model::set_space_dim - wrong space dimension." );

                    // set space dimension
                    mSpaceDim = aSpaceDim;
                }

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags specific to certain material models
                 */
                virtual void reset_specific_eval_flags(){};

                //------------------------------------------------------------------------------
                /**
                 * set material model dof types
                 * @param[ in ] aDofTypes a list of group of dof types
                 */
                void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * set material model dof types
                 * @param[ in ] aDofTypes a list of group of dof types
                 * @param[ in ] aDofStrings a list of strings to describe the dof types
                 */
                virtual void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                        moris::Cell< std::string >                  aDofStrings )
                {
                    MORIS_ERROR( false, "Material_Model::set_dof_type_list - Not implemented for base class." );
                }

                //------------------------------------------------------------------------------
                /**
                 * return a cell of dof types
                 * @param[ in ] mDofTypes a cell of cell of dof types
                 */
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list() const
                {
                    return mDofTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * build a map for the dof types
                 */
                void build_dof_type_map();

                //------------------------------------------------------------------------------
                /**
                 * get dof type map
                 * @param[ out ] mDofTypeMap map for the dof types
                 */
                const Matrix< DDSMat > & get_dof_type_map()
                {
                    return mDofTypeMap;
                }

                //------------------------------------------------------------------------------
                /**
                 * set field interpolator manager
                 * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
                 */
                void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager );

                //------------------------------------------------------------------------------
                /**
                 * get properties
                 * @param[ out ] mProperties cell of property pointers
                 */
                moris::Cell< std::shared_ptr< Property > > & get_properties()
                {
                    return mProperties;
                };

                //------------------------------------------------------------------------------
                /**
                 * create a global dof type list including material and property dependencies
                 */
                void build_global_dof_type_list();

                //------------------------------------------------------------------------------
                /**
                 * initialize storage variables and evaluation flags specific to some child MMs
                 * function is called in the build_global_dof_type_list()
                 */
                virtual void initialize_spec_storage_vars_and_eval_flags(){};

                //------------------------------------------------------------------------------
                /**
                 * get global dof type list
                 * @param[ out ] mGlobalDofTypes global list of dof type
                 */
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list();

                //------------------------------------------------------------------------------
                /**
                 * get non unique list of dof type for the material model
                 * @param[ in ] aDofTypes a cell of dof type to fill
                 */
                void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * build global dof type map
                 */
                void build_global_dof_type_map();

                //------------------------------------------------------------------------------
                /**
                 * get global dof type map
                 * @param[ out ] mGlobalDofTypeMap global map for the dof types
                 */
                const Matrix< DDSMat > & get_global_dof_type_map();

                //------------------------------------------------------------------------------
                /**
                 * check dependency on a given group of dof types
                 * @param[ in ]  aDofType       a group of dof types
                 * @param[ out ] tDofDependency a bool true if dependency on dof type
                 *
                 */
                bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType );

                //------------------------------------------------------------------------------
                // SPECIFIC INTERNAL ENERGY (FIRST EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the specific internal energy and its derivatives
                 */
                virtual void eval_Eint()
                {
                    MORIS_ERROR( false, " Material_Model::eval_Eint - This function does nothing. " );
                };

                virtual void eval_EintDot()
                {
                    MORIS_ERROR( false, " Material_Model::eval_EintDot - This function does nothing. " );
                };

                virtual void eval_dEintdx()
                {
                    MORIS_ERROR( false, " Material_Model::eval_dnEintdxn - This function does nothing. " );
                };

                virtual void eval_d2Eintdx2()
                {
                    MORIS_ERROR( false, " Material_Model::eval_dnEintdxn - This function does nothing. " );
                };               

                //------------------------------------------------------------------------------
                /**
                 * evaluate the specific internal energy derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_EintDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_EintDOF - This function does nothing. " );
                };

                virtual void eval_EintDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_EintDotDOF - This function does nothing. " );
                };

                virtual void eval_dEintdxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes)
                {
                    MORIS_ERROR( false, " Material_Model::eval_dEintdxDOF - This function does nothing. " );
                };

                virtual void eval_d2Eintdx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes)
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Eintdx2DOF - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the specific internal energy and its derivatives
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnEintdxn specific internal energy and derivatives
                 */
                const Matrix< DDRMat > & Eint();
                const Matrix< DDRMat > & EintDot();
                const Matrix< DDRMat > & dnEintdxn( uint aOrder );

                //------------------------------------------------------------------------------
                /**
                 * get the specific internal energy and its derivatives wrt to the DoF types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnEintdxnDOF specific internal energy and derivatives
                 */
                const Matrix< DDRMat > & EintDOF( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & EintDotDOF( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnEintdxnDOF( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );

                //------------------------------------------------------------------------------
                // DENSITY (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the density and its derivatives
                 */
                virtual void eval_density()
                {
                    MORIS_ERROR( false, " Material_Model::eval_density - This function does nothing. " );
                };

                virtual void eval_DensityDot()
                {
                    MORIS_ERROR( false, " Material_Model::eval_DensityDot - This function does nothing. " );
                };

                virtual void eval_dDensitydx()
                {
                    MORIS_ERROR( false, " Material_Model::eval_dDensitydx - This function does nothing. " );
                };

                virtual void eval_d2Densitydx2()
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Densitydx2 - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic density derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_DensityDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_DensityDOF - This function does nothing. " );
                };

                virtual void eval_DensityDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_DensityDotDOF - This function does nothing. " );
                };

                virtual void eval_dDensitydxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_dDensitydxDOF - This function does nothing. " );
                };

                virtual void eval_d2Densitydx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Densitydx2DOF - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic density and its derivatives
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnDensitydxn thermodynamic density and derivatives
                 */
                const Matrix< DDRMat > & density(){ return ( this->*m_get_Density )(); };
                const Matrix< DDRMat > & DensityDot(){ return ( this->*m_get_DensityDot )(); };
                const Matrix< DDRMat > & dnDensitydxn( uint aOrder ){ return ( this->*m_get_dnDensitydxn )( aOrder ); };

                const Matrix< DDRMat > & density_dep();
                const Matrix< DDRMat > & DensityDot_dep();
                const Matrix< DDRMat > & dnDensitydxn_dep( uint aOrder );      
                
                const Matrix< DDRMat > & density_triv();
                const Matrix< DDRMat > & DensityDot_triv();
                const Matrix< DDRMat > & dnDensitydxn_triv( uint aOrder );                          

                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic density and its derivatives wrt to the DoF types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnDensitydxnDOF thermodynamic density and derivatives
                 */
                const Matrix< DDRMat > & DensityDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_DensityDof )( aDofType );
                };
                const Matrix< DDRMat > & DensityDotDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_DensityDotDof )( aDofType ); 
                };
                const Matrix< DDRMat > & dnDensitydxnDOF( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder )
                { 
                    return ( this->*m_get_dnDensitydxnDof )( aDofType, aOrder ); 
                };        

                const Matrix< DDRMat > & DensityDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & DensityDotDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnDensitydxnDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );     

                const Matrix< DDRMat > & DensityDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & DensityDotDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnDensitydxnDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );                                            

                //------------------------------------------------------------------------------
                // PRESSURE (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the pressure and its derivatives
                 */
                virtual void eval_pressure()
                {
                    MORIS_ERROR( false, " Material_Model::eval_pressure - This function does nothing. " );
                };

                virtual void eval_PressureDot()
                {
                    MORIS_ERROR( false, " Material_Model::eval_PressureDot - This function does nothing. " );
                };

                virtual void eval_dPressuredx()
                {
                    MORIS_ERROR( false, " Material_Model::eval_dPressuredx - This function does nothing. " );
                };

                virtual void eval_d2Pressuredx2()
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Pressuredx2 - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic pressure derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_PressureDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_PressureDOF - This function does nothing. " );
                };

                virtual void eval_PressureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_PressureDotDOF - This function does nothing. " );
                };

                virtual void eval_dPressuredxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_dPressuredxDOF - This function does nothing. " );
                };

                virtual void eval_d2Pressuredx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Pressuredx2DOF - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic pressure and its derivatives
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnPressuredxn thermodynamic pressure and derivatives
                 */
                const Matrix< DDRMat > & pressure(){ return ( this->*m_get_Pressure )(); };
                const Matrix< DDRMat > & PressureDot(){ return ( this->*m_get_PressureDot )(); };
                const Matrix< DDRMat > & dnPressuredxn( uint aOrder ){ return ( this->*m_get_dnPressuredxn )( aOrder ); };

                const Matrix< DDRMat > & pressure_dep();
                const Matrix< DDRMat > & PressureDot_dep();
                const Matrix< DDRMat > & dnPressuredxn_dep( uint aOrder );      
                
                const Matrix< DDRMat > & pressure_triv();
                const Matrix< DDRMat > & PressureDot_triv();
                const Matrix< DDRMat > & dnPressuredxn_triv( uint aOrder );                          

                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic pressure and its derivatives wrt to the DoF types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnPressuredxnDOF thermodynamic pressure and derivatives
                 */
                const Matrix< DDRMat > & PressureDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_PressureDof )( aDofType );
                };
                const Matrix< DDRMat > & PressureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_PressureDotDof )( aDofType ); 
                };
                const Matrix< DDRMat > & dnPressuredxnDOF( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder )
                { 
                    return ( this->*m_get_dnPressuredxnDof )( aDofType, aOrder ); 
                };        

                const Matrix< DDRMat > & PressureDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & PressureDotDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnPressuredxnDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );     

                const Matrix< DDRMat > & PressureDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & PressureDotDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnPressuredxnDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );                

                //------------------------------------------------------------------------------
                // TEMPERATURE (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the temperature and its derivatives
                 */
                virtual void eval_temperature()
                {
                    MORIS_ERROR( false, " Material_Model::eval_temperature - This function does nothing. " );
                };

                virtual void eval_TemperatureDot()
                {
                    MORIS_ERROR( false, " Material_Model::eval_TemperatureDot - This function does nothing. " );
                };

                virtual void eval_dTemperaturedx()
                {
                    MORIS_ERROR( false, " Material_Model::eval_dTemperaturedx - This function does nothing. " );
                };

                virtual void eval_d2Temperaturedx2()
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Temperaturedx2 - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic temperature derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_TemperatureDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_TemperatureDOF - This function does nothing. " );
                };

                virtual void eval_TemperatureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_TemperatureDotDOF - This function does nothing. " );
                };

                virtual void eval_dTemperaturedxDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_dTemperaturedxDOF - This function does nothing. " );
                }; 

                virtual void eval_d2Temperaturedx2DOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_d2Temperaturedx2DOF - This function does nothing. " );
                };               

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic temperature and its derivatives
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnTemperaturedxn thermodynamic temperature and derivatives
                 */
                const Matrix< DDRMat > & temperature(){ return ( this->*m_get_Temperature )(); };
                const Matrix< DDRMat > & TemperatureDot(){ return ( this->*m_get_TemperatureDot )(); };
                const Matrix< DDRMat > & dnTemperaturedxn( uint aOrder ){ return ( this->*m_get_dnTemperaturedxn )( aOrder ); };

                const Matrix< DDRMat > & temperature_dep();
                const Matrix< DDRMat > & TemperatureDot_dep();
                const Matrix< DDRMat > & dnTemperaturedxn_dep( uint aOrder );      
                
                const Matrix< DDRMat > & temperature_triv();
                const Matrix< DDRMat > & TemperatureDot_triv();
                const Matrix< DDRMat > & dnTemperaturedxn_triv( uint aOrder );                          

                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic temperature and its derivatives wrt to the DoF types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] mdnTemperaturedxnDOF thermodynamic temperature and derivatives
                 */
                const Matrix< DDRMat > & TemperatureDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_TemperatureDof )( aDofType ); 
                };
                const Matrix< DDRMat > & TemperatureDotDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
                { 
                    return ( this->*m_get_TemperatureDotDof )( aDofType ); 
                };
                const Matrix< DDRMat > & dnTemperaturedxnDOF( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder )
                { 
                    return ( this->*m_get_dnTemperaturedxnDof )( aDofType, aOrder ); 
                };        

                const Matrix< DDRMat > & TemperatureDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & TemperatureDotDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnTemperaturedxnDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder );     

                const Matrix< DDRMat > & TemperatureDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & TemperatureDotDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType );
                const Matrix< DDRMat > & dnTemperaturedxnDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder ); 

                //------------------------------------------------------------------------------
                // THERMODYNAMIC QUANTITIES
                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic quantities
                 */
                virtual void eval_VolumeExpansivity()
                {
                    MORIS_ERROR( false, " Material_Model::eval_VolumeExpansivity - This function does nothing. " );
                };

                virtual void eval_IsothermalCompressibility()
                {
                    MORIS_ERROR( false, " Material_Model::eval_IsothermalCompressibility - This function does nothing. " );
                };

                virtual void eval_Cv()
                {
                    MORIS_ERROR( false, " Material_Model::eval_Cv - This function does nothing. " );
                };

                virtual void eval_Cp()
                {
                    MORIS_ERROR( false, " Material_Model::eval_Cp - This function does nothing. " );
                };   

                virtual void eval_Gamma()
                {
                    MORIS_ERROR( false, " Material_Model::eval_Gamma - This function does nothing. " );
                };
            
                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic quantity derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_VolumeExpansivityDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_VolumeExpansivityDOF - This function does nothing. " );
                };

                virtual void eval_IsothermalCompressibilityDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_IsothermalCompressibilityDOF - This function does nothing. " );
                };

                virtual void eval_CvDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_CvDOF - This function does nothing. " );
                };

                virtual void eval_CpDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_CpDOF - This function does nothing. " );
                };   

                virtual void eval_GammaDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Material_Model::eval_GammaDOF - This function does nothing. " );
                };

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic quantities
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] tQuantity thermodynamic quantity
                 */
                const Matrix< DDRMat > & AlphaP();
                const Matrix< DDRMat > & BetaT();

                const Matrix< DDRMat > & Cv(){ return ( this->*m_get_Cv )(); };
                const Matrix< DDRMat > & Cp(){ return ( this->*m_get_Cp )(); };
                const Matrix< DDRMat > & Gamma(){ return ( this->*m_get_Gamma )(); };

                const Matrix< DDRMat > & Cv_dep();
                const Matrix< DDRMat > & Cp_dep();
                const Matrix< DDRMat > & Gamma_dep();     
  
                const Matrix< DDRMat > & Cv_triv();
                const Matrix< DDRMat > & Cp_triv();
                const Matrix< DDRMat > & Gamma_triv();   

                //------------------------------------------------------------------------------
                /**
                 * get the thermodynamic quantities and their derivatives wrt to the DoF types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the derivative (only for x-derivs)
                 * @param[ out ] aQuantityDOF thermodynamic quantity and derivatives
                 */
                const Matrix< DDRMat > & AlphaPDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                const Matrix< DDRMat > & BetaTDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                const Matrix< DDRMat > & CvDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                { 
                    return ( this->*m_get_CvDof )( aDofTypes ); 
                };
                const Matrix< DDRMat > & CpDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                { 
                    return ( this->*m_get_CpDof )( aDofTypes ); 
                };
                const Matrix< DDRMat > & GammaDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                { 
                    return ( this->*m_get_GammaDof )( aDofTypes ); 
                };

                const Matrix< DDRMat > & CvDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                const Matrix< DDRMat > & CpDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                const Matrix< DDRMat > & GammaDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes );     
                
                const Matrix< DDRMat > & CvDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                const Matrix< DDRMat > & CpDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                const Matrix< DDRMat > & GammaDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes );                 

                //------------------------------------------------------------------------------
                // FINITE DIFFERENCE FUNCTIONS FOR COMPUTING DOF DERIVATIVES
                //------------------------------------------------------------------------------
                /**
                 * evaluate the specific internal energy wrt to a dof type
                 * @param[ in ] aDofTypes       dof type wrt which the derivative is evaluated
                 * @param[ in ] aEintDOF_FD     matrix to fill with derivative of the specific internal energy
                 * @param[ in ] aPerturbation   real to perturb for FD
                 * @param[ in ] aDerivOrder     order of the spatial derivatives ( includes 0 )
                 * @param[ in ] aFDSchemeType   enum for FD scheme
                 */
                void eval_EintDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aEintDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5 );

                void eval_EintDotDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aEintDotDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5 );

                void eval_dnEintdxnDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adnEintdxnDOF_FD,
                        real                                 aPerturbation,
                        uint                                 aOrder,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5 );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the third dependent thermodynamic variable wrt to a dof type
                 * @param[ in ] aDofTypes       dof type wrt which the derivative is evaluated
                 * @param[ in ] aTDvarDOF_FD    matrix to fill with derivative of the 3rd thermodynamic variable
                 * @param[ in ] aPerturbation   real to perturb for FD
                 * @param[ in ] aDerivOrder     order of the spatial derivatives ( includes 0 )
                 * @param[ in ] aFDSchemeType   enum for FD scheme
                 */
                void eval_TDvarDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aTDvarDOF_FD,
                        real                                 aPerturbation,
                        MSI::Dof_Type                        aTDvar,
                        fem::FDScheme_Type                   aFDSchemeType );

                void eval_TDvarDotDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aTDvarDotDOF_FD,
                        real                                 aPerturbation,
                        MSI::Dof_Type                        aTDvar,
                        fem::FDScheme_Type                   aFDSchemeType );

                void eval_dnTDvardxnDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adnTDvardxnDOF_FD,
                        real                                 aPerturbation,
                        MSI::Dof_Type                        aTDvar,
                        uint                                 aOrder,
                        fem::FDScheme_Type                   aFDSchemeType );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic quantities wrt to a dof type
                 * @param[ in ] aDofTypes        dof type wrt which the derivative is evaluated
                 * @param[ in ] aQuantityDOF_FD  matrix to fill with derivative of the specific internal energy
                 * @param[ in ] aQuantityString  string specifying which quantity should be finite differenced
                 * @param[ in ] aPerturbation    real to perturb for FD
                 * @param[ in ] aFDSchemeType    enum for FD scheme
                 */
                void eval_QuantityDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aQuantityDOF_FD,
                        std::string                          aQuantityString,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5 );                 

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_MATERIAL_MODEL_HPP_ */
