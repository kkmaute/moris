/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Dof_Manager.hpp
 *
 */

#ifndef SRC_FEM_CL_DOF_MANAGER_HPP_
#define SRC_FEM_CL_DOF_MANAGER_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace MSI
    {
        class Pdof_Host;
        class Model_Solver_Interface;
        class Equation_Set;
        class Dof_Manager
        {
            private:
                moris::Cell < Pdof_Host * >          mPdofHostList;  // List of all pdof hosts
                moris::Cell < Adof * >               mAdofList;      // List of all adofs
                moris::Cell<moris::Cell < Adof * > > mAdofListOwned; // List of all owned adofs
                moris::Cell<moris::Cell < Adof * > > mAdofListOwnedAndShared; // List of all owned adofs
                //moris::Cell < Adof * >       mAdofListShared;

                moris::Cell< enum Dof_Type >   mPdofTypeList;          // List containing all used unique dof types.
                Matrix< DDSMat >               mPdofTypeMap;           // Map which maps the unique dof types onto consecutive values.

                Matrix< DDUMat >               mPdofHostTimeLevelList; // List containing the number of time levels per dof type.
                Matrix< DDUMat >               mTimeLevelOffsets;      // Type offsets  created through time levels
                Matrix< IdMat >                mCommTable;             // Communication table. As and input from the model.
                Model_Solver_Interface       * mModelSolverInterface;  // Model solver interface pointer

                Cell< moris::map< moris::moris_id, moris::moris_index > > mAdofGlobaltoLocalMap;
                moris::sint mNumMaxAdofs = -1;
                moris::uint mNumOwnedAdofs = 0;

                Matrix< DDSMat > mTypeTimeIndentifierToTypeMap;

                bool mUseHMR = false;

                // List containing the number of time levels per dof type.
                // FIXME
                Matrix< DDUMat > mTimePerDofType;

                //------------------------------------------------------------------------------
                /**
                 * @brief Returns the maximal number of pdof hosts. This function is tested by the test [Dof_Manager_Max_Pdof_Host]
                 *
                 * @param[in] aListEqnObj   List containing all the equation objects.
                 *
                 */
                moris::uint initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj );

                //------------------------------------------------------------------------------
                /**
                 * @brief Communicates the local list with pdof types. This function is called inside initialize_pdof_type_list.
                 *
                 * @param[in] aPdofTypeList   List of pdof types
                 *
                 */
                void communicate_dof_types( moris::Cell< enum Dof_Type > & aPdofTypeList );

                //------------------------------------------------------------------------------
                /**
                 * @brief Creates a map which maps the dof type enum values to a consecutive list of dof type indices. These indices are used in the rest of the module.
                 * This function is called inside initialize_pdof_type_list.
                 *
                 * @param[in] aPdofTypeList   List of pdof types
                 *
                 */
                void create_dof_type_map();

                //------------------------------------------------------------------------------
                /**
                 * @brief Creates a list containing the number of time levels per pdof type. This function is tested by the test [Dof_Manager_Pdof_Host_Time_Level]
                 *
                 */
                void initialize_pdof_host_time_level_list();

                //------------------------------------------------------------------------------
                /**
                 * @brief Communicates the local time levels over all processors. Should only have an effect if a dof type is not existent on a processor.
                 * This function is called by initialize_pdof_host_time_level_list.
                 *
                 * @param[in] aTimeLevelList   List of number of time levels per dof type.
                 *
                 */
                void communicate_time_list( Matrix< DDUMat > & aTimeLevelList );

                //------------------------------------------------------------------------------
                /**
                 * @brief This member function checks if a adof which is shared has a corresponding owned adof on the owning processor.
                 * If the owning processor never created the adof, it will create a dummy adof and flag it as owned.
                 * This function is called inside create_adofs.
                 *
                 * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
                 *
                 */
                void
                communicate_check_if_owned_adof_exists(
                        moris::Cell< moris::Cell < Adof * > > & tAdofListofTypes,
                        Matrix< IndexMat > const              & aDiscretizationIndexPerTypeAndTime);

                //------------------------------------------------------------------------------
                /**
                 * @brief This member function calculates and communicated the adof Id offset for each processor. This function is called inside create_adofs.
                 *
                 * @param[in] aNumOwnedAdofs   Number of owned adofs
                 *
                 */
                moris::uint communicate_adof_offsets( const moris::uint & aNumOwnedAdofs );

                //------------------------------------------------------------------------------
                /**
                 * @brief In this function communicates the onwed adof Ids to the shared adofs.. This function is called inside create_adofs.
                 *
                 * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
                 * @param[in] aListSharedAdofIds A list containing the shared adof Ids.
                 * @param[in] aListSharedAdofPos A list containing the adof position of the shared adof Id
                 *
                 */
                void communicate_shared_adof_ids(
                        moris::Cell< moris::Cell < Adof * > > const & aAdofListofTypes,
                        Matrix< IndexMat >                    const & aDiscretizationIndexPerTypeAndTime,
                        moris::Cell< Matrix< DDUMat > >             & aListSharedAdofIds,
                        moris::Cell< Matrix< DDUMat > >             & aListSharedAdofPos);

                //------------------------------------------------------------------------------
                /**
                 * @brief This functon determines the maximal adof index
                 *
                 * @param[in] aMaxAdofInd Maximal possible adof index
                 *
                 */
                void get_max_adof_ind( moris::sint & aMaxAdofInd );

                //------------------------------------------------------------------------------
                /**
                 * @brief Function calling the routine to set the owned Adof Ids. Options are ordered by type or ordered by host
                 *
                 * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
                 * @param[in] aAdofOffsets     Adof offsets for this processor.
                 *
                 */
                void set_owned_adofs_ids( const moris::Cell<moris::Cell < Adof * > > & aAdofListofTypes,
                        const uint                                 & aAdofOffsets );

                //------------------------------------------------------------------------------
                /**
                 * @brief Function setting the Adof Ids per host
                 *
                 * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
                 * @param[in] aAdofOffsets     Adof offsets for this processor.
                 *
                 */
                void set_owned_adofs_ids_by_type( const moris::Cell<moris::Cell < Adof * > > & aAdofListofTypes,
                        const uint                                 & aAdofOffsets );

                //------------------------------------------------------------------------------
                /**
                 * @brief Function setting the Adof Ids per host/underlying interpolation basis. Note this only is useful if all adofs have the same interpolation order.
                 *        If multiple interpolation meshes are used the underlying basis might not be the same although they might share the same index.
                 *
                 * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
                 * @param[in] aAdofOffsets     Adof offsets for this processor.
                 *
                 */
                void set_owned_adofs_ids_by_host( const moris::Cell<moris::Cell < Adof * > > & aAdofListofTypes,
                        const uint                                 & aAdofOffsets );

                //------------------------------------------------------------------------------

                void get_descretization_index_for_adof_list_of_types(
                        Matrix< IndexMat > & tDiscretizationIndexPerTypeAndTime );

            public:
                Dof_Manager(){};

                //-----------------------------------------------------------------------------------------------------------
                Dof_Manager( const Matrix< IdMat > aCommTable ) : mCommTable( aCommTable )
                {
                    mUseHMR = true;
                };

                //-----------------------------------------------------------------------------------------------------------
                Dof_Manager( const Matrix< IdMat > aCommTable,
                        Model_Solver_Interface * aModelSolverInterface ) : mCommTable( aCommTable ),
                                mModelSolverInterface( aModelSolverInterface )
                {
                    mUseHMR = true;
                };

                //-----------------------------------------------------------------------------------------------------------
                ~Dof_Manager();

                //-----------------------------------------------------------------------------------------------------------
                void set_adof_map( const moris::map< moris::moris_id, moris::moris_index > * aAdofLocaltoGlobalMap )
                {
                    mAdofGlobaltoLocalMap.resize( 1 );
                    mAdofGlobaltoLocalMap( 0 ) = *aAdofLocaltoGlobalMap;
                }

                //-----------------------------------------------------------------------------------------------------------
                Cell< moris::map< moris::moris_id, moris::moris_index > > & get_adof_map(  )
                    {
                    return mAdofGlobaltoLocalMap;
                    }

                //-----------------------------------------------------------------------------------------------------------

                void set_max_num_adofs( const moris::sint & aNumMaxAdofs )
                {
                    mNumMaxAdofs = aNumMaxAdofs;
                }

                //-----------------------------------------------------------------------------------------------------------
                /**
                 * @brief Initializes a list with all used physical dof types. This function is tested by the test [Dof_Mgn_create_unique_dof_type_list][Dof_Mgn_create_unique_dof_type_map_matrix]
                 *
                 * @param[in] aListEqnObj   List containing all the equation objects.
                 *
                 */
                void initialize_pdof_type_list( Cell< MSI::Equation_Set * > & aListEqnBlock );

                //-----------------------------------------------------------------------------------------------------------

                void set_time_levels_for_type( const enum Dof_Type aDofType,
                        const moris::uint   aNumTimeLevels )
                {
                    moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aDofType ) );

                    MORIS_ASSERT( tDofTypeIndex != -1, "Dof_Manager::set_time_levels_for_type(), Requested dof type does not exist");

                    mTimePerDofType( tDofTypeIndex ) = aNumTimeLevels;
                };

                //-----------------------------------------------------------------------------------------------------------

                moris::uint get_time_levels_for_type( const enum Dof_Type aDofType )
                {
                    moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aDofType ) );

                    return mTimePerDofType( tDofTypeIndex );
                };

                //-----------------------------------------------------------------------------------------------------------

                const Matrix< DDUMat > & get_time_levels( )
                    {
                    return mTimePerDofType;
                    };

                //-----------------------------------------------------------------------------------------------------------
                /**
                 * @brief Initializes list with pdof hosts. This function is tested by the test   [Dof_Mgn_ini_pdof_host_list]
                 *
                 * @param[in] aListEqnObj   List containing all the equation objects.
                 *
                 */
                void initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj );

                //-----------------------------------------------------------------------------------------------------------
                /**
                 * @brief This member function creates the list containing all active adofs. This function is tested by the test [ Dof_Mgn_create_adofs ]
                 * [Dof_Mgn_create_adofs_parallell_1][Dof_Mgn_create_adofs_parallell_2]
                 *
                 */
                void create_adofs();

                //-----------------------------------------------------------------------------------------------------------
                /**
                 * @brief This member function tells each pdof to int the pdof host list to get the T-matrix entries. This function is tested by the test [Dof_Mgn_Set_T_Matrix]
                 *
                 */
                void set_pdof_t_matrix();

                //-----------------------------------------------------------------------------------------------------------
                /**
                 * @brief This member function returns the size of the adof list.
                 *
                 *@param[out] Number of adofs
                 *
                 */
                moris::uint get_num_owned_adofs()
                {
                    return mNumOwnedAdofs;
                };

                //-----------------------------------------------------------------------------------------------------------

                moris::sint get_pdof_index_for_type( enum Dof_Type aDofType)
                {
                    moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aDofType ) );

                    MORIS_ASSERT( tDofTypeIndex != -1,"get_pdof_index_for_type() requested dof does not exist in map");

                    return tDofTypeIndex;
                };

                //-----------------------------------------------------------------------------------------------------------

                moris::Matrix< DDSMat > get_unique_adof_mesh_indices();

                const moris::Matrix< DDSMat > & get_typetime_identifier_to_type_map()
                    {
                    return mTypeTimeIndentifierToTypeMap;
                    };

                //-----------------------------------------------------------------------------------------------------------

                Cell< moris::Cell < Adof * > > get_owned_adofs()
                    {
                    return mAdofListOwned;
                    };

                //-----------------------------------------------------------------------------------------------------------
                Matrix< DDSMat > get_local_adof_ids();

                //-----------------------------------------------------------------------------------------------------------

                Matrix< DDSMat > get_local_adof_ids( const moris::Cell< enum Dof_Type > & aListOfDofTypes );

                //-----------------------------------------------------------------------------------------------------------

                Matrix< DDSMat > get_local_overlapping_adof_ids( const moris::Cell< enum Dof_Type > & aListOfDofTypes );

                //-----------------------------------------------------------------------------------------------------------
                Matrix< DDSMat > get_local_overlapping_adof_ids();

                //-----------------------------------------------------------------------------------------------------------

                //this function is for HMR use only. It creates a map between MSI adof inds and HMR adof inds
                Matrix< DDUMat > get_adof_ind_map();

                //-----------------------------------------------------------------------------------------------------------

                enum Dof_Type get_dof_type_enum( moris::uint aDofType )
                {
                    return mPdofTypeList( aDofType );
                };

                //-----------------------------------------------------------------------------------------------------------

                moris::uint get_num_dof_types()
                {
                    MORIS_ASSERT( mPdofTypeList.size() != 0, "Dof_Manager::get_num_dof_types(), pdof type list not initialized yet");

                    return mPdofTypeList.size();
                };

        };
    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_DOF_MANAGER_HPP_ */

