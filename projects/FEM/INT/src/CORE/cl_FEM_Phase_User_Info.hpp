/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Phase_User_Info.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_
#define SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_

//MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
namespace moris
{
    namespace fem
    {

        class Constitutive_Model;
        class Stabilization_Parameter;

        //------------------------------------------------------------------------------
        /**
         * Phase_User_Info
         */
        class Phase_User_Info
        {
            protected :

                // phase name
                std::string mPhaseName;

                // phase index
                moris::Matrix< moris::IndexMat > mPhaseIndex;

                // phase dof type list
                moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

                // matrix to check if a dof type was already set to the phase
                Matrix< DDSMat > mDofCheck;

                // phase dv type list
                moris::Cell< moris::Cell< PDV_Type > > mDvTypes;

                // matrix to check if a pdv type was already set to the phase
                Matrix< DDSMat > mPdvCheck;

                // constitutive models
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mCMs;
                uint mNumCMs = 0;

                // constitutive model map
                std::map< std::string, uint > mCMMap;

                // material models
                moris::Cell< std::shared_ptr< fem::Material_Model > > mMMs;
                uint mNumMMs = 0;

                // constitutive model map
                std::map< std::string, uint > mMMMap;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * trivial constructor
                 */
                Phase_User_Info()
                {
                    // set size for mDofCheck
                    mDofCheck.set_size( static_cast< uint >( MSI::Dof_Type::END_ENUM ), 1, -1 );

                    // set size for mPdvCheck
                    mPdvCheck.set_size( static_cast< uint >( PDV_Type::UNDEFINED ), 1, -1 );
                }

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~Phase_User_Info(){};

                //------------------------------------------------------------------------------
                /**
                 * print names
                 */
                void print_names()
                {
                    // print the phase name
                    std::cout<<"Phase name: "<<mPhaseName<<std::endl;

                    // print the phase indices
                    for( uint iPhaseIndex = 0; iPhaseIndex < mPhaseIndex.numel(); iPhaseIndex++ )
                    {
                        std::cout << "Phase index: " << mPhaseIndex( iPhaseIndex ) << std::endl;
                    }

                    // print CM names
                    for ( uint iCM = 0; iCM < mCMs.size(); iCM++ )
                    {
                        std::cout << "CM name: " << mCMs( iCM )->get_name() << std::endl;
                    }

                    // print mM names
                    for ( uint iMM = 0; iMM < mMMs.size(); iMM++ )
                    {
                        std::cout << "MM name: " << mMMs( iMM )->get_name() << std::endl;
                    }

                }

                //------------------------------------------------------------------------------
                /**
                 * set the phase name
                 * @param[ in ] aPhaseName mesh set name
                 */
                void set_phase_name( std::string aPhaseName )
                {
                    mPhaseName = aPhaseName;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the phase index
                 * @param[ out ] mPhaseName mesh set name
                 */
                std::string get_phase_name()
                {
                    return mPhaseName;
                }

                //------------------------------------------------------------------------------
                /**
                 * set phase mesh indices
                 * @param[ in ] aPhaseIndex cell of mesh indices for phase
                 */
                void set_phase_indices( moris::Matrix< moris::IndexMat > aPhaseIndex )
                {
                    mPhaseIndex = aPhaseIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * get phase mesh indices
                 * @param[ out ] mPhaseIndex cell of mesh indices for phase
                 */
                moris::Matrix< moris::IndexMat > & get_phase_indices()
                {
                    return mPhaseIndex;
                }

                //------------------------------------------------------------------------------
                /**
                 * set dof type list
                 * @param[ in ] aDofTypes list of group of dof types
                 */
                void set_dof_type_list(
                        const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
                {
                    mDofTypes = aDofTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * add dof type to list
                 * @param[ in ] aDofTypes group of dof types to add to list of dof type
                 */
                void add_dof_type_to_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
                {
                    // loop over all dof types
                    for ( uint iType = 0; iType < aDofTypes.size(); ++iType )
                    {
                        // get dof type index in enum list
                        uint tDofIndex = static_cast< uint >( aDofTypes( iType )( 0 ) );

                        // if dof type not added to phase dof type list
                        if( mDofCheck( tDofIndex ) == -1 )
                        {
                            // add dof type group to dof type list
                            mDofTypes.push_back( aDofTypes( iType ) );

                            // set check to 1
                            mDofCheck( tDofIndex ) = 1;
                        }
                    }
                }

                //------------------------------------------------------------------------------

                const Matrix< DDSMat > & get_dof_type_check_list()
                {
                    return mDofCheck;
                }

                //------------------------------------------------------------------------------
                /**
                 * get dof type list
                 * @param[ out ] mDofTypes list of group of dof types
                 */
                const
                moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list()
                {
                    return mDofTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * set dv type list
                 * @param[ out ] mDvTypes list of group of dv types
                 */
                void set_dv_type_list( const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes )
                {
                    mDvTypes = aDvTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * get dv type list
                 * @param[ out ] mDvTypes list of group of dv types
                 */
                const
                moris::Cell< moris::Cell< PDV_Type > > & get_dv_type_list()
                {
                    return mDvTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * add pdv type to list
                 * @param[ in ] aPdvTypes group of pdv types to add to list of pdv type
                 */
                void add_pdv_type_to_list( moris::Cell< PDV_Type > & aPdvTypes )
                {
                    // get pdv type index in enum list
                    uint tPdvIndex = static_cast< uint >( aPdvTypes( 0 ) );

                    // if pdv type not added to phase pdv type list
                    if( mPdvCheck( tPdvIndex ) == -1 )
                    {
                        // add pdv type group to pdv type list
                        mDvTypes.push_back( aPdvTypes );

                        // set check to 1
                        mPdvCheck( tPdvIndex ) = 1;
                    }
                }

                Matrix< DDSMat > & get_pdv_type_check_list()
                {
                    return mPdvCheck;
                }

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * set CMs
                 * @param[ in ] aCMs list of CM pointers
                 */
                void set_CMs( const moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & aCMs )
                {
                    mCMs = aCMs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set CM
                 * @param[ in ] aCM constitutive model pointer
                 */
                void set_CM( std::shared_ptr< fem::Constitutive_Model > aCM )
                {
                    // get CM name
                    std::string tCMName = aCM->get_name();

                    // if CM was not set before
                    if ( mCMMap.find( tCMName ) == mCMMap.end() )
                    {
                        // add CM to list
                        mCMs.push_back( aCM );

                        // add CM index to map
                        mCMMap[ tCMName ] = mNumCMs++;
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * get CM by name
                 * @param[ out ] aCMName constitutive model name requested
                 */
                std::shared_ptr< fem::Constitutive_Model > get_CM_by_name( std::string aCMName )
                {
                    // check for unknown CM name
                    MORIS_ERROR( mCMMap.find( aCMName ) != mCMMap.end(),
                            "Phase_User_Info::get_CM_by_name - Unknown aCMName for %s : %s \n",
                            mPhaseName.c_str(),
                            aCMName.c_str() );

                    // return CM for name
                    return mCMs( mCMMap[ aCMName ] );
                }

                //------------------------------------------------------------------------------
                /**
                 * get CMs
                 * @param[ out ] mCMs list of constitutive model pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & get_CMs() const
                {
                    return mCMs;
                }

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * set MMs
                 * @param[ in ] aMMs list of MM pointers
                 */
                void set_MMs( const moris::Cell< std::shared_ptr< fem::Material_Model > > & aMMs )
                {
                    mMMs = aMMs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set MM
                 * @param[ in ] aMM material model pointer
                 */
                void set_MM( std::shared_ptr< fem::Material_Model > aMM )
                {
                    // get MM name
                    std::string tMMName = aMM->get_name();

                    // if MM was not set before
                    if ( mMMMap.find( tMMName ) == mMMMap.end() )
                    {
                        // add MM to list
                        mMMs.push_back( aMM );

                        // add MM index to map
                        mMMMap[ tMMName ] = mNumMMs++;
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * get MM by name
                 * @param[ out ] aMMName material model name requested
                 */
                std::shared_ptr< fem::Material_Model > get_MM_by_name( std::string aMMName )
                {
                    // check for unknown MM name
                    MORIS_ERROR( mMMMap.find( aMMName ) != mMMMap.end(),
                            "Phase_User_Info::get_MM_by_name - Unknown aMMName for %s : %s \n",
                            mPhaseName.c_str(),
                            aMMName.c_str() );

                    // return MM for name
                    return mMMs( mMMMap[ aMMName ] );
                }

                //------------------------------------------------------------------------------
                /**
                 * get MMs
                 * @param[ out ] mMMs list of material model pointers
                 */
                const
                moris::Cell< std::shared_ptr< fem::Material_Model > > & get_MMs() const
                {
                    return mMMs;
                }

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PHASE_USER_INFO_HPP_ */

