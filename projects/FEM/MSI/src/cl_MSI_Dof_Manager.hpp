/*
 * cl_Dof_Manager.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
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
//------------------------------------------------------------------------------
        class Pdof_Host;
        //class Adof;
//------------------------------------------------------------------------------

        class Dof_Manager
        {
//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
            moris::Cell < Pdof_Host * >  mPdofHostList;           // List of all pdof hosts
            moris::Cell < Adof * >       mAdofList;               // List of all adofs
            moris::Cell < Adof * >       mAdofListOwned;          // List of all owned adofs
            //moris::Cell < Adof * >       mAdofListShared;

            moris::Cell< enum Dof_Type > mPdofTypeList;            // List containing all used unique dof types.
            Matrix< DDSMat >    mPdofTypeMap;             // Map which maps the unique dof types onto consecutive values.
            Matrix< DDUMat >    mPdofHostTimeLevelList;   // List containing the number of time levels per dof type.
            Matrix< IdMat >     mCommTable;               // Communication table. As and input from the model.

            const moris::map< moris::moris_id, moris::moris_index >  * mAdofGlobaltoLocalMap = nullptr;
            moris::sint mNumMaxAdofs = -1;

            bool mUseHMR = false;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Dof_Manager(){};

//-----------------------------------------------------------------------------------------------------------

            Dof_Manager( const Matrix< IdMat > aCommTable ) : mCommTable( aCommTable )
            {
                mUseHMR = true;
            };

//-----------------------------------------------------------------------------------------------------------

            ~Dof_Manager();

//-----------------------------------------------------------------------------------------------------------

            void
            set_adof_map( const moris::map< moris::moris_id, moris::moris_index > * aAdofLocaltoGlobalMap )
            {
                mAdofGlobaltoLocalMap = aAdofLocaltoGlobalMap;
            }

//-----------------------------------------------------------------------------------------------------------

            void
            set_max_num_adofs( const moris::sint & aNumMaxAdofs )
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
            void
            initialize_pdof_type_list( moris::Cell < Equation_Object* > & aListEqnObj );

//-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Initializes list with pdof hosts. This function is tested by the test   [Dof_Mgn_ini_pdof_host_list]
             *
             * @param[in] aListEqnObj   List containing all the equation objects.
             *
             */
            void
            initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj );

//-----------------------------------------------------------------------------------------------------------

            /**
             * @brief This member function creates the list containing all active adofs. This function is tested by the test [ Dof_Mgn_create_adofs ]
             * [Dof_Mgn_create_adofs_parallell_1][Dof_Mgn_create_adofs_parallell_2]
             *
             */
            void
            create_adofs();

//-----------------------------------------------------------------------------------------------------------
            /**
             * @brief This member function tells each pdof to int the pdof host list to get the T-matrix entries. This function is tested by the test [Dof_Mgn_Set_T_Matrix]
             *
             */
            void
            set_pdof_t_matrix();
//-----------------------------------------------------------------------------------------------------------

            /**
             * @brief This member function returns the size of the adof list.
             *
             *@param[out] Number of adofs
             *
             */
            moris::uint get_num_adofs() { return mAdofListOwned.size(); };

//-----------------------------------------------------------------------------------------------------------

            Matrix< DDSMat >
            get_local_adof_ids();

//-----------------------------------------------------------------------------------------------------------

            Matrix< DDSMat >
            get_local_overlapping_adof_ids();

//-----------------------------------------------------------------------------------------------------------

            //this function is for HMR use only. It creates a map between MSI adof inds and HMR adof inds
            Matrix< DDUMat >
            get_adof_ind_map();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * @brief Returns the maximal number of pdof hosts. This function is tested by the test [Dof_Manager_Max_Pdof_Host]
             *
             * @param[in] aListEqnObj   List containing all the equation objects.
             *
             */
            moris::uint
            initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj );

//------------------------------------------------------------------------------

            /**
             * @brief Communicates the local list with pdof types. This function is called inside initialize_pdof_type_list.
             *
             * @param[in] aPdofTypeList   List of pdof types
             *
             */
            void
            communicate_dof_types( moris::Cell< enum Dof_Type > & aPdofTypeList );

//------------------------------------------------------------------------------

            /**
             * @brief Creates a map which maps the dof type enum values to a consecutive list of dof type indices. These indices are used in the rest of the module.
             * This function is called inside initialize_pdof_type_list.
             *
             * @param[in] aPdofTypeList   List of pdof types
             *
             */
            void
            create_dof_type_map();

//------------------------------------------------------------------------------

            /**
             * @brief Creates a list containing the number of time levels per pdof type. This function is tested by the test [Dof_Manager_Pdof_Host_Time_Level]
             *
             */
            void
            initialize_pdof_host_time_level_list();

//------------------------------------------------------------------------------

            /**
             * @brief Communicates the local time levels over all processors. Should only have an effect if a dof type is not existent on a processor.
             * This function is called by initialize_pdof_host_time_level_list.
             *
             * @param[in] aTimeLevelList   List of number of time levels per dof type.
             *
             */
            void
            communicate_time_list( Matrix< DDUMat > & aTimeLevelList );

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
            communicate_check_if_owned_adof_exists( moris::Cell< moris::Cell < Adof * > > & tAdofListofTypes );

//------------------------------------------------------------------------------

            /**
             * @brief This member function calculates and communicated the adof Id offset for each processor. This function is called inside create_adofs.
             *
             * @param[in] aNumOwnedAdofs   Number of owned adofs
             *
             */
            moris::uint
            communicate_adof_offsets( const moris::uint & aNumOwnedAdofs );

//------------------------------------------------------------------------------

            /**
             * @brief In this function communicates the onwed adof Ids to the shared adofs.. This function is called inside create_adofs.
             *
             * @param[in] tAdofListofTypes A temporary list containing a list of adofs for every doftype and timelevel.
             * @param[in] aListSharedAdofIds A list containing the shared adof Ids.
             * @param[in] aListSharedAdofPos A list containing the adof position of the shared adof Id
             *
             */
            void communicate_shared_adof_ids(const moris::Cell< moris::Cell < Adof * > > & aAdofListofTypes,
                    Matrix< DDUMat >             & aListSharedAdofIds,
                    Matrix< DDUMat >             & aListSharedAdofPos);
        };
//------------------------------------------------------------------------------

    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_DOF_MANAGER_HPP_ */
