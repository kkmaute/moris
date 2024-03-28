/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Dof_Manager.cpp
 *
 */

#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Pdof_Host.hpp"

#include "fn_sort.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace MSI
    {
        Dof_Manager::~Dof_Manager()
        {
            for ( uint Ik = 0; Ik < mPdofHostList.size(); Ik++ )
            {
                delete mPdofHostList( Ik );
            }
            mPdofHostList.clear();
            for ( uint Ik = 0; Ik < mAdofList.size(); Ik++ )
            {
                delete mAdofList( Ik );
            }
            mAdofList.clear();
        }

        //-----------------------------------------------------------------------------------------------------------
        uint
        Dof_Manager::initialize_max_number_of_possible_pdof_hosts( Vector< Equation_Object* >& aListEqnObj )
        {
            // Ask how many equation objects
            uint tNumEqnObj       = aListEqnObj.size();
            uint tMaxNumPdofHosts = 0;

            // loop over all equation objects, asking for their number of pdof hosts
            for ( uint Ii = 0; Ii < tNumEqnObj; Ii++ )
            {
                tMaxNumPdofHosts = std::max( tMaxNumPdofHosts, aListEqnObj( Ii )->get_max_pdof_hosts_ind() );
            }

            MORIS_ERROR( tMaxNumPdofHosts < MORIS_SINT_MAX, "Invalid Index in PDOF host" );

            // Add 1 because c++ is 0 based
            tMaxNumPdofHosts = tMaxNumPdofHosts + 1;

            return tMaxNumPdofHosts;
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Dof_Manager::initialize_pdof_type_list( Vector< MSI::Equation_Set* >& aListEqnBlock )
        {
            // Reserve of temporary pdof type list
            Vector< enum Dof_Type > tTemporaryPdofTypeList;

            tTemporaryPdofTypeList.reserve( static_cast< int >( Dof_Type::END_ENUM ) + 1 );

            Matrix< DDUMat > tListToCheckIfEnumExist( ( static_cast< int >( Dof_Type::END_ENUM ) + 1 ), 1, 0 );

            // Get number of equation objects
            uint tNumEquationBlocks = aListEqnBlock.size();

            // loop over all equation objects, asking for their pdof types
            for ( uint Ii = 0; Ii < tNumEquationBlocks; Ii++ )
            {
                // Ask equation object for its dof types
                const Vector< enum Dof_Type >& tDofType = aListEqnBlock( Ii )->get_unique_dof_type_list();

                // Loop over all dof types
                for ( uint Ik = 0; Ik < tDofType.size(); Ik++ )
                {
                    // Set 1 at position of the enum value
                    tListToCheckIfEnumExist( static_cast< int >( tDofType( Ik ) ) ) = 1;
                }
            }

            // Loop over tListToCheckIfEnumExist. Add Ddof type to list if corresponding value = 1
            for ( uint Ij = 0; Ij < tListToCheckIfEnumExist.numel(); Ij++ )
            {
                if ( tListToCheckIfEnumExist( Ij ) == 1 )
                {
                    tTemporaryPdofTypeList.push_back( (Dof_Type)Ij );
                }
            }

            // Shrink pdof type list to fit
            tTemporaryPdofTypeList.shrink_to_fit();

            // Communicate dof types so that all processors have the same unique list
            this->communicate_dof_types( tTemporaryPdofTypeList );

            // Create a map
            this->create_dof_type_map();
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::communicate_dof_types( Vector< enum Dof_Type >& aPdofTypeList )
        {
            // Get processor size
            int tSize = par_size();

            // Get number of local dof types
            sint tNumLocalDofTypes = aPdofTypeList.size();

            // Get number of global dof types
            sint tNumMaxGlobalDofTypes = sum_all( tNumLocalDofTypes );

            if ( par_rank() == 0 )
            {
                // Set size of of pdof type list = number of global types
                mPdofTypeList.resize( tNumMaxGlobalDofTypes );
            }

            // Create list containing the number of local dof types
            Vector< sint > tNumLocalDofTypesList( tSize );

            // Insert number of local dof types into list containing the number of local dof types
            MPI_Allgather( &tNumLocalDofTypes, 1, MPI_UNSIGNED, ( tNumLocalDofTypesList.data() ).data(), 1, MPI_UNSIGNED, MPI_COMM_WORLD );

            // Create list containing the offsets of the local dof types in relation to processor 0
            Vector< sint > tDofTypeOffset( tSize );

            // Fill the list with the corresponding offsets
            for ( int Ip = 1; Ip < tSize; ++Ip )
            {
                tDofTypeOffset( Ip ) = tDofTypeOffset( Ip - 1 ) + tNumLocalDofTypesList( Ip - 1 );
            }

            // Assemble list containing all used dof types. Dof types are not unique
            MPI_Gatherv(
                    ( ( aPdofTypeList.data() ).data() ),
                    tNumLocalDofTypes,
                    MPI_UNSIGNED,
                    ( mPdofTypeList.data() ).data(),
                    ( tNumLocalDofTypesList.data() ).data(),
                    ( tDofTypeOffset.data() ).data(),
                    MPI_UNSIGNED,
                    0,
                    MPI_COMM_WORLD );

            // Temporary variable for mPdofTypeList size
            uint tPdofTypeListSize;

            if ( par_rank() == 0 )
            {
                // Sort this created list
                std::sort( ( mPdofTypeList.data() ).data(), ( mPdofTypeList.data() ).data() + mPdofTypeList.size() );

                // use std::unique and std::distance to create list containing all used dof types. This list is unique
                auto last = std::unique( ( mPdofTypeList.data() ).data(), ( mPdofTypeList.data() ).data() + mPdofTypeList.size() );
                auto pos  = std::distance( ( mPdofTypeList.data() ).data(), last );

                mPdofTypeList.resize( pos );

                tPdofTypeListSize = mPdofTypeList.size();
            }

            // Bcast size of mPdofTypeList on processor 0
            MPI_Bcast( &tPdofTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

            // Resize mPdofTypeList on all processors
            mPdofTypeList.resize( tPdofTypeListSize );

            // Bcast unique mPdofTypeList to all processors
            MPI_Bcast( ( mPdofTypeList.data() ).data(), mPdofTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );

            mTimePerDofType.set_size( tPdofTypeListSize, 1, 1 );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::create_dof_type_map()
        {
            // Get number of unique adofs of this equation object
            uint tNumUniquePdofTypes = mPdofTypeList.size();

            // Get maximal dof type enum number
            sint tMaxDofTypeEnumNumber = 0;

            // Loop over all pdof types to get the highest enum index
            for ( uint Ii = 0; Ii < tNumUniquePdofTypes; Ii++ )
            {
                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mPdofTypeList( Ii ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;

            // Set size of maping matrix
            mPdofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

            // Loop over all pdof types to create the mapping matrix
            for ( uint Ii = 0; Ii < tNumUniquePdofTypes; Ii++ )
            {
                mPdofTypeMap( static_cast< int >( mPdofTypeList( Ii ) ) ) = Ii;
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::initialize_pdof_host_list( Vector< Equation_Object* >& aListEqnObj )
        {
            uint tNumUsedDofTypes = mPdofTypeList.size();
            uint tNumEquationObj  = aListEqnObj.size();
            uint tMaxNumPdofHosts = this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );

            // reserve the maximal needed amount of memory for a temporary PdofHostList
            Vector< Pdof_Host* > tPdofHostList( tMaxNumPdofHosts, nullptr );

            // Loop over all equation objects. Ask them to create their pdof hosts. Pdof hosts are stored in tPdofHostList
            for ( uint Ii = 0; Ii < tNumEquationObj; Ii++ )
            {
                aListEqnObj( Ii )->create_my_pdof_hosts( tNumUsedDofTypes, mPdofTypeMap, mTimePerDofType, tPdofHostList );
            }

            // Determine number of Pdof Hosts
            uint tNumPdofHosts = 0;
            for ( uint Ik = 0; Ik < tMaxNumPdofHosts; Ik++ )
            {
                // If pointer in temporary pdof host list exists add one to number of pdof hosts
                if ( tPdofHostList( Ik ) != NULL )
                {
                    tNumPdofHosts = tNumPdofHosts + 1;
                }
            }

            // in case of repeated update of the equation sets, we have to delete the old pdof hosts
            for ( uint Ik = 0; Ik < mPdofHostList.size(); Ik++ )
            {
                delete mPdofHostList( Ik );
            }
            // Set size of List containing all pdof hosts
            mPdofHostList.resize( tNumPdofHosts );

            uint counter = 0;
            // add pointers to pdof hosts into list of pdof hosts
            for ( uint Ij = 0; Ij < tMaxNumPdofHosts; Ij++ )
            {
                // If pointer in temporary pdof host list exists add one to number of pdof hosts
                if ( tPdofHostList( Ij ) != NULL )
                {
                    mPdofHostList( counter ) = tPdofHostList( Ij );
                    counter                  = counter + 1;
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::initialize_pdof_host_time_level_list()
        {
            uint tNumPdofTypes = mPdofTypeList.size();
            uint tNumPdofHosts = mPdofHostList.size();

            // Set size of the list of pdof host time levels
            Matrix< DDUMat > tPdofHostTimeLevelList( tNumPdofTypes, 1, 0 );

            // Loop over pdof types
            for ( uint Ik = 0; Ik < tNumPdofTypes; Ik++ )
            {
                // loop over pdof hosts
                for ( uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
                {
                    uint tNumTimeLevels = mPdofHostList( Ii )->get_num_time_levels_of_type( Ik );

                    if ( tNumTimeLevels != 0 )
                    {
                        // Get the number of time levels per pdof type
                        tPdofHostTimeLevelList( Ik ) = tNumTimeLevels;

                        // Break loop if number of time levels for this pdof type is set.
                        break;
                    }
                }
            }

            // communicate time level and set mPdofHostTimeLevelList
            this->communicate_time_list( tPdofHostTimeLevelList );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::communicate_time_list( Matrix< DDUMat >& aTimeLevelList )
        {
            mPdofHostTimeLevelList.set_size( mPdofTypeList.size(), 1, 0 );

            // Get the maximal value for all time list values
            for ( uint Ii = 0; Ii < aTimeLevelList.numel(); Ii++ )
            {
                mPdofHostTimeLevelList( Ii ) = max_all( aTimeLevelList( Ii ) );
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::communicate_check_if_owned_adof_exists(
                Vector< Vector< Adof* > >& aAdofListofTypes,
                Matrix< IndexMat > const & aDiscretizationIndexPerTypeAndTime )
        {
            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap( mCommTable.max() + 1, 1, -1 );

            uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ) ) = Ik;
            }

            // Loop over all different adof types and time in this temporary list
            for ( uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
            {
                Vector< Matrix< DDUMat > > tSharedAdofPosGlobal( tNumCommProcs );

                // Set Mat to store number of shared adofs per processor
                Matrix< DDUMat > tNumSharedAdofsPerProc( tNumCommProcs, 1, 0 );

                // Loop over adofs per type/time. Count number of adofs per proc which have to be communicated
                for ( uint Ib = 0; Ib < aAdofListofTypes( Ij ).size(); Ib++ )
                {
                    // Check if adof at this position is not NULL
                    if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                    {
                        // Check if owning processor is this processor
                        if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                        {
                            // get owning processor
                            moris::moris_id tProcID = aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                            sint tProcIdPos = tCommTableMap( tProcID );

                            MORIS_ASSERT( tProcIdPos != -1, "Dof_Manager::communicate_check_if_owned_adof_exists: Map returns proc rank -1. Check communication table" );

                            // Add +1 to the processor number of shared dofs per processor
                            tNumSharedAdofsPerProc( tProcIdPos )++;
                        }
                    }
                }

                // Set size of the moris::Mats in the Cell
                for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    if ( tNumSharedAdofsPerProc( Ik ) != 0 )
                    {
                        tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik ), 1 );
                    }
                }

                // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
                Matrix< DDUMat > tShredAdofPosPerProc( tNumCommProcs, 1, 0 );

                // Loop over adofs per type
                for ( uint Ia = 0; Ia < aAdofListofTypes( Ij ).size(); Ia++ )
                {
                    // Check if adof at this position is not NULL
                    if ( aAdofListofTypes( Ij )( Ia ) != NULL )
                    {
                        // Check if owning processor is this processor
                        if ( aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                        {
                            // Get owning processor
                            uint tProcID = aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor();

                            sint tProcIdPos = tCommTableMap( tProcID );

                            // Add owning processor id to moris::Mat
                            tSharedAdofPosGlobal( tProcIdPos )( tShredAdofPosPerProc( tProcIdPos ), 0 ) =
                                    aAdofListofTypes( Ij )( Ia )->get_adof_external_id();

                            tShredAdofPosPerProc( tProcIdPos ) = tShredAdofPosPerProc( tProcIdPos ) + 1;
                        }
                    }
                }

                // receiving list
                Vector< Matrix< DDUMat > > tMatsToReceive;

                barrier();

                // Communicate position of shared adofs to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedAdofPosGlobal,
                        tMatsToReceive );

                // Loop over all Mats set dummy owned adofs
                for ( uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    for ( uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        moris_index tDiscretizationIndex = aDiscretizationIndexPerTypeAndTime( Ij );

                        // Get owned adof Id
                        uint tLocalAdofInd = mAdofGlobaltoLocalMap( tDiscretizationIndex ).find( tMatsToReceive( Ik )( Ii ) );

                        if ( aAdofListofTypes( Ij )( tLocalAdofInd ) == NULL )
                        {
                            aAdofListofTypes( Ij )( tLocalAdofInd ) = new Adof();
                            aAdofListofTypes( Ij )( tLocalAdofInd )->set_adof_owning_processor( par_rank() );

                            // Set external adof index. Used for HMR ordering
                            aAdofListofTypes( Ij )( tLocalAdofInd )->set_adof_external_ind( tLocalAdofInd );
                        }
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        uint
        Dof_Manager::communicate_adof_offsets( const uint& aNumOwnedAdofs )
        {
            MORIS_LOG_SPEC( "Total number of DOFs", sum_all( aNumOwnedAdofs ) );

            // Get list containing the number of owned adofs of each processor
            Matrix< DDUMat > tNumOwnedAdofsList;
            comm_gather_and_broadcast( aNumOwnedAdofs, tNumOwnedAdofsList );

            Matrix< DDUMat > tOwnedAdofsOffsetList( tNumOwnedAdofsList.numel(), 1, 0 );

            // Loop over all entries to create the offsets. Starting with 1
            for ( uint Ij = 1; Ij < tOwnedAdofsOffsetList.numel(); Ij++ )
            {
                // Add the number of owned adofs of the previous processor to the offset of the previous processor
                tOwnedAdofsOffsetList( Ij, 0 ) = tOwnedAdofsOffsetList( Ij - 1, 0 ) + tNumOwnedAdofsList( Ij - 1, 0 );
            }

            return tOwnedAdofsOffsetList( par_rank(), 0 );
        }

        //----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::communicate_shared_adof_ids(
                Vector< Vector< Adof* > > const & aAdofListofTypes,
                Matrix< IndexMat > const &        aDiscretizationIndexPerTypeAndTime,
                Vector< Matrix< DDUMat > >&       aListSharedAdofIds,
                Vector< Matrix< DDUMat > >&       aListSharedAdofPos )
        {
            // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
            Matrix< DDSMat > tCommTableMap( mCommTable.max() + 1, 1, -1 );

            uint tNumCommProcs = mCommTable.numel();

            // Loop over communication table to fill the communication table map
            for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                tCommTableMap( mCommTable( Ik ) ) = Ik;
            }

            // Loop over all different adof types and time in this temporary list
            for ( uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
            {
                uint tCounter       = 0;
                uint tSharedCounter = 0;

                Vector< Matrix< DDUMat > > tSharedAdofPosGlobal( tNumCommProcs );
                Vector< Matrix< DDUMat > > tSharedAdofPosLocal( tNumCommProcs );

                // Set Mat to store number of shared adofs per processor
                Matrix< DDUMat > tNumSharedAdofsPerProc( tNumCommProcs, 1, 0 );

                // Loop over adofs per type/time
                for ( uint Ib = 0; Ib < aAdofListofTypes( Ij ).size(); Ib++ )
                {
                    // Check if adof at this position is not NULL
                    if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                    {
                        // Check if owning processor is this processor
                        if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                        {
                            // get owning procssor
                            moris::moris_id tProcID = aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                            sint tProcIdPos = tCommTableMap( tProcID );

                            // Add +1 to the processor number of shared dofs per processor
                            tNumSharedAdofsPerProc( tProcIdPos )++;

                            tSharedCounter++;
                        }
                    }
                }

                // Set size of the moris::Mats in the Cell
                for ( uint Ik = 0; Ik < tNumCommProcs; Ik++ )
                {
                    if ( tNumSharedAdofsPerProc( Ik ) != 0 )
                    {
                        tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik ), 1 );
                        tSharedAdofPosLocal( Ik ).set_size( tNumSharedAdofsPerProc( Ik ), 1 );
                    }
                }

                // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
                Matrix< DDUMat > tShredAdofPosPerProc( tNumCommProcs, 1, 0 );

                // Loop over adofs per type
                for ( uint Ia = 0; Ia < aAdofListofTypes( Ij ).size(); Ia++ )
                {
                    // Check if adof at this position is not NULL
                    if ( aAdofListofTypes( Ij )( Ia ) != NULL )
                    {
                        // Check if owning processor is this processor
                        if ( aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                        {
                            // Get owning processor
                            moris::moris_id tProcID = aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor();

                            sint tProcIdPos = tCommTableMap( tProcID );

                            // Add owning processor id to moris::Mat
                            tSharedAdofPosGlobal( tProcIdPos )( tShredAdofPosPerProc( tProcIdPos ) ) =
                                    aAdofListofTypes( Ij )( Ia )->get_adof_external_id();

                            // Add adof position to Mat
                            tSharedAdofPosLocal( tProcIdPos )( tShredAdofPosPerProc( tProcIdPos ) ) = tCounter;

                            tShredAdofPosPerProc( tProcIdPos )++;
                        }
                    }
                    tCounter++;
                }

                // receiving list
                Vector< Matrix< DDUMat > > tMatsToReceive;

                barrier();

                // Communicate position of shared adofs to the owning processor
                communicate_mats(
                        mCommTable,
                        tSharedAdofPosGlobal,
                        tMatsToReceive );

                // Create List of Mats containing the shared node Ids
                Vector< Matrix< DDUMat > > tSharesAdofIdList( tNumCommProcs );

                // Loop over all Mats setting the size
                for ( uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    tSharesAdofIdList( Ik ).set_size( tMatsToReceive( Ik ).numel(), 1 );
                }

                // Loop over all received positions and get the adof id of the owning adof
                for ( uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
                {
                    for ( uint Ii = 0; Ii < tMatsToReceive( Ik ).numel(); Ii++ )
                    {
                        moris_index tDiscretizationIndex = aDiscretizationIndexPerTypeAndTime( Ij );

                        // Get owned adof Id
                        uint tLocalAdofInd = mAdofGlobaltoLocalMap( tDiscretizationIndex ).find( tMatsToReceive( Ik )( Ii ) );

                        MORIS_ASSERT( ( aAdofListofTypes( Ij )( tLocalAdofInd )->get_adof_owning_processor() ) == par_rank(),
                                "Dof_Manager::communicate_shared_adof_ids: Adof not owned by this processor" );

                        tSharesAdofIdList( Ik )( Ii ) = aAdofListofTypes( Ij )( tLocalAdofInd )->get_adof_id();

                        // tSharesAdofIdList( Ik )( Ii, 0 ) = ( aAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) )->get_adof_id();
                    }
                }

                Vector< Matrix< DDUMat > > tMatsToReceive2;

                barrier();

                // Communicate owned adof Id back to the processor with the shared adof
                communicate_mats( mCommTable,
                        tSharesAdofIdList,
                        tMatsToReceive2 );

                uint tAdofPosCounter = 0;

                aListSharedAdofIds( Ij ).set_size( tSharedCounter, 1, MORIS_UINT_MAX );
                aListSharedAdofPos( Ij ).set_size( tSharedCounter, 1, MORIS_UINT_MAX );

                // assemble Ids in list of shared adof ids and assemble the corresponding postions
                for ( uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
                {
                    if ( tMatsToReceive2( Ik ).numel() >= 1 )
                    {
                        aListSharedAdofIds( Ij )( { tAdofPosCounter, tAdofPosCounter + tMatsToReceive2( Ik ).numel() - 1 }, { 0, 0 } ) =
                                tMatsToReceive2( Ik ).matrix_data();

                        aListSharedAdofPos( Ij )( { tAdofPosCounter, tAdofPosCounter + tSharedAdofPosLocal( Ik ).numel() - 1 }, { 0, 0 } ) =
                                tSharedAdofPosLocal( Ik ).matrix_data();

                        tAdofPosCounter = tAdofPosCounter + tMatsToReceive2( Ik ).numel();
                    }
                }

                if ( aListSharedAdofIds( Ij ).numel() != 0 )
                {
                    MORIS_ASSERT( aListSharedAdofIds( Ij ).max() != MORIS_UINT_MAX,
                            "Dof_Manager::communicate_shared_adof_ids(), communicated Ids not set correctly" );

                    MORIS_ASSERT( aListSharedAdofPos( Ij ).max() != MORIS_UINT_MAX,
                            "Dof_Manager::communicate_shared_adof_ids(), positions for communicated Ids not set correctly" );
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::get_max_adof_ind( sint& aMaxAdofInd )
        {
            // Get max entry of node adof if pdof host list exists
            if ( mNumMaxAdofs == -1 )
            {
                uint tNumPdofHosts = mPdofHostList.size();

                sint tMaxNumAdofMeshInd = mModelSolverInterface->get_max_adof_index();

                moris::Matrix< DDSMat > tAdofMeshIndexExists( tMaxNumAdofMeshInd, 1, -1 );
                moris::Matrix< DDSMat > tAdofMeshIndex( tMaxNumAdofMeshInd, 1, -1 );
                uint                    tIndexCounter = 0;

                // Loop over all dof types. Determine the adof orders of these dof types and put them into a list. can be 1,2,3
                for ( uint Ij = 0; Ij < mPdofTypeList.size(); Ij++ )
                {
                    // Ask for adof order for this dof type
                    uint tAdofMeshInd = (uint)mModelSolverInterface->get_adof_index_for_type( Ij );

                    if ( tAdofMeshIndexExists( tAdofMeshInd ) == -1 )
                    {
                        tAdofMeshIndexExists( tAdofMeshInd ) = 1;

                        tAdofMeshIndex( tIndexCounter++ ) = tAdofMeshInd;
                    }
                }
                tAdofMeshIndex.resize( tIndexCounter, 1 );

                if ( tNumPdofHosts != 0 )
                {
                    for ( uint Ia = 0; Ia < tAdofMeshIndex.numel(); Ia++ )
                    {
                        uint tAdofMeshInd = tAdofMeshIndex( Ia );

                        for ( uint Ik = 0; Ik < tNumPdofHosts; Ik++ )
                        {
                            moris::fem::Node_Base* tNode = mPdofHostList( Ik )->get_node_obj_ptr();

                            aMaxAdofInd = std::max( aMaxAdofInd, ( tNode->get_adof_indices( tAdofMeshInd ) ).max() );
                        }
                    }
                }
                // Add one because c++ is 0 based. ==> List size has to be tMaxNodeAdofId + 1
                aMaxAdofInd = aMaxAdofInd + 1;
            }
            else if ( mNumMaxAdofs >= 0 )
            {
                aMaxAdofInd = mNumMaxAdofs;
            }
            else
            {
                MORIS_ERROR( false, "MSI::Dof_Manager: Check number of adofs" );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Dof_Manager::set_owned_adofs_ids(
                const Vector< Vector< Adof* > >& aAdofListofTypes,
                const uint&                      aAdofOffsets )
        {
            Parameter_List& tParameterlist = mModelSolverInterface->get_msi_parameterlist();

            if ( tParameterlist.get< bool >( "order_adofs_by_host" ) )
            {
                MORIS_LOG_INFO( "Setting adofs by host" );
                this->set_owned_adofs_ids_by_host( aAdofListofTypes, aAdofOffsets );
            }
            else
            {
                this->set_owned_adofs_ids_by_type( aAdofListofTypes, aAdofOffsets );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Dof_Manager::set_owned_adofs_ids_by_type(
                const Vector< Vector< Adof* > >& aAdofListofTypes,
                const uint&                      aAdofOffsets )
        {
            uint tCounterId = aAdofOffsets;
            // uint tCounterShared = 0;
            uint tCounter = 0;

            // loop over temporary adof list. Add pointers to adofs into list of adofs
            for ( uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
            {
                uint tCounterOwned = 0;

                uint tCounterOwnedAndShared = 0;

                for ( uint Ib = 0; Ib < aAdofListofTypes( Ij ).size(); Ib++ )
                {
                    // If pointer in temporary adofs list exists, add adof
                    if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                    {
                        if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() == par_rank() )
                        {
                            // Add adof to owned adof list
                            mAdofListOwned( Ij )( tCounterOwned ) = aAdofListofTypes( Ij )( Ib );

                            // Set adof Id
                            mAdofListOwned( Ij )( tCounterOwned )->set_adof_id( tCounterId );

                            tCounterOwned++;
                            tCounterId++;
                        }

                        // Add adof to owned  and shared adof list
                        mAdofListOwnedAndShared( Ij )( tCounterOwnedAndShared++ ) = aAdofListofTypes( Ij )( Ib );

                        // Add adof to adof list
                        mAdofList( tCounter++ ) = aAdofListofTypes( Ij )( Ib );
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Dof_Manager::set_owned_adofs_ids_by_host(
                const Vector< Vector< Adof* > >& aAdofListofTypes,
                const uint&                      aAdofOffsets )
        {
            // FIXME by host only makes sense if all adofs have the same interpolation index
            // otherwise it will be just mixed weirdly

            uint tCounterId = aAdofOffsets;
            // uint tCounterShared = 0;
            uint tCounter = 0;

            // Counter for every woned type and time
            moris::Matrix< DDUMat > tCounterOwned( aAdofListofTypes.size(), 1, 0 );
            // Counter for every owned and shared type and time
            moris::Matrix< DDUMat > tCounterOwnedAndShared( aAdofListofTypes.size(), 1, 0 );

            // loop over temporary adof list. Add pointers to adofs into list of adofs
            for ( uint Ib = 0; Ib < aAdofListofTypes( 0 ).size(); Ib++ )
            {
                for ( uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
                {
                    // If pointer in temporary adofs list exists, add adof
                    if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                    {
                        if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() == par_rank() )
                        {
                            // Add adof to owned adof list
                            mAdofListOwned( Ij )( tCounterOwned( Ij ) ) = aAdofListofTypes( Ij )( Ib );

                            // Set adof Id
                            mAdofListOwned( Ij )( tCounterOwned( Ij ) )->set_adof_id( tCounterId );

                            tCounterOwned( Ij )++;
                            tCounterId++;
                        }

                        // Add adof to owned and shared adof list
                        mAdofListOwnedAndShared( Ij )( tCounterOwnedAndShared( Ij )++ ) = aAdofListofTypes( Ij )( Ib );

                        // Add adof to adof list
                        mAdofList( tCounter++ ) = aAdofListofTypes( Ij )( Ib );
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::get_descretization_index_for_adof_list_of_types(
                Matrix< IndexMat >& tDiscretizationIndexPerTypeAndTime )
        {
            uint tCounter = 0;

            // loop over all dof types
            for ( uint Ik = 0; Ik < mPdofTypeList.size(); Ik++ )
            {
                // get discretization index for this dof type
                uint tDiscretizationIndex =
                        (uint)mModelSolverInterface->get_adof_index_for_type( Ik );

                // loop over all time levels for this dof type
                for ( uint Ii = 0; Ii < mPdofHostTimeLevelList( Ik ); Ii++ )
                {
                    // assign discretization index  to all time levels for this dof type
                    tDiscretizationIndexPerTypeAndTime( tCounter++ ) = tDiscretizationIndex;
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Dof_Manager::create_adofs()
        {
            // set time levels and populate mPdofHostTimeLevelList
            this->initialize_pdof_host_time_level_list();

            mTimeLevelOffsets.set_size( mPdofHostTimeLevelList.numel(), 1, 0 );
            for ( uint Ik = 1; Ik < mPdofHostTimeLevelList.numel(); Ik++ )
            {
                mTimeLevelOffsets( Ik ) = mTimeLevelOffsets( Ik - 1 ) + mPdofHostTimeLevelList( Ik - 1 );
            }

            // Get number of pdof types and size of pdof host list
            uint tNumPdofHosts      = mPdofHostList.size();
            uint tNumTypeTimeLevels = sum( mPdofHostTimeLevelList );
            sint tMaxNodeAdofId     = -1;

            this->get_max_adof_ind( tMaxNodeAdofId );

            // Create temporary Vector containing lists of temporary adofs
            Vector< Vector< Adof* > > tAdofListofTypes( tNumTypeTimeLevels );

            // build discretization index list for types and time.
            // Every entry of this list corresponds to a adof list in tAdofListofTypes
            Matrix< IndexMat > tDiscretizationIndexPerTypeAndTime( tNumTypeTimeLevels, 1 );

            // fill list with correct discretization indices for type and time
            this->get_descretization_index_for_adof_list_of_types( tDiscretizationIndexPerTypeAndTime );

            for ( uint Ik = 0; Ik < tNumTypeTimeLevels; Ik++ )
            {
                tAdofListofTypes( Ik ).resize( tMaxNodeAdofId, nullptr );
            }

            // Loop over all pdof hosts and get the adofs
            for ( uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
            {
                mPdofHostList( Ii )->get_adofs( mTimeLevelOffsets, tAdofListofTypes, mModelSolverInterface, mUseHMR );
            }

            // Check if shared adof exists
            if ( par_size() > 1 )
            {
                this->communicate_check_if_owned_adof_exists( tAdofListofTypes, tDiscretizationIndexPerTypeAndTime );
            }

            // Determine number of owned and shared adofs
            uint tNumAdofs = 0;

            // Multigrid type and time identifier;
            sint tAdofTypeTimeIdentifier = 0;
            mTypeTimeIndentifierToTypeMap.set_size( tAdofListofTypes.size(), 1, -1 );

            moris::Matrix< DDUMat > tNumOwnedAdofsPerTypeTime( tNumTypeTimeLevels, 1, 0 );

            moris::Matrix< DDUMat > tNumOwnedAndShardAdofsPerTypeTime( tNumTypeTimeLevels, 1, 0 );

            uint tCounterTypeTime = 1;

            // Loop over all adofs determine the total number and the number of owned ones
            for ( uint Ik = 0; Ik < tAdofListofTypes.size(); Ik++ )
            {
                for ( uint Ia = 0; Ia < tAdofListofTypes( Ik ).size(); Ia++ )
                {
                    // If pointer in temporary adof list exists. Add one to number of owned adofs
                    if ( ( tAdofListofTypes( Ik )( Ia ) != NULL ) )
                    {
                        tNumAdofs++;

                        // If owning processor equals this processor then its an owned adof
                        if ( ( tAdofListofTypes( Ik )( Ia )->get_adof_owning_processor() == par_rank() ) )
                        {
                            tNumOwnedAdofsPerTypeTime( Ik )++;
                        }

                        tNumOwnedAndShardAdofsPerTypeTime( Ik )++;

                        // Add type/time identifier to Adof
                        tAdofListofTypes( Ik )( Ia )->set_adof_type_time_identifier( tAdofTypeTimeIdentifier );
                    }
                }

                // Multigrid Type time identifier to type map  -> Type for Type time identifier
                if ( tCounterTypeTime < mTimeLevelOffsets.numel() )
                {
                    if ( Ik < mTimeLevelOffsets( tCounterTypeTime ) )
                    {
                        mTypeTimeIndentifierToTypeMap( Ik ) = tCounterTypeTime - 1;
                    }
                    else if ( Ik == mTimeLevelOffsets( tCounterTypeTime ) )
                    {
                        tCounterTypeTime++;
                        mTypeTimeIndentifierToTypeMap( Ik ) = tCounterTypeTime - 1;
                    }
                    else
                    {
                        MORIS_ASSERT( false, "Dof_Manager::create_adofs(), " );
                    }
                }
                else
                {
                    mTypeTimeIndentifierToTypeMap( Ik ) = tCounterTypeTime - 1;
                }

                tAdofTypeTimeIdentifier++;
            }

            // Calculate number of owned and shared adofs
            mNumOwnedAdofs = sum( tNumOwnedAdofsPerTypeTime );
            //        uint tNumSharedAdofs = tNumAdofs - mNumOwnedAdofs;

            // Get adof offset for this processor
            uint tAdofOffset = this->communicate_adof_offsets( mNumOwnedAdofs );

            // remove previously allocated adofs (if any)
            for ( uint Ik = 0; Ik < mAdofList.size(); Ik++ )
            {
                delete mAdofList( Ik );
            }
            // Set size of List containing all adofs
            mAdofList.resize( tNumAdofs, nullptr );
            mAdofListOwned.resize( tNumTypeTimeLevels );
            mAdofListOwnedAndShared.resize( tNumTypeTimeLevels );

            for ( uint Ik = 0; Ik < mAdofListOwned.size(); Ik++ )
            {
                mAdofListOwned( Ik ).resize( tNumOwnedAdofsPerTypeTime( Ik ), nullptr );
                mAdofListOwnedAndShared( Ik ).resize( tNumOwnedAndShardAdofsPerTypeTime( Ik ), nullptr );
            }
            // mAdofListShared.resize( tNumSharedAdofs );

            // set the owned adof ids
            this->set_owned_adofs_ids( tAdofListofTypes, tAdofOffset );

            if ( par_size() > 1 )
            {
                // Create vector storing information about adof ids and position for communication
                Vector< Matrix< DDUMat > > tListSharedAdofIds( tNumTypeTimeLevels );
                Vector< Matrix< DDUMat > > tListSharedAdofPos( tNumTypeTimeLevels );

                // Communicate shared-owned adof ids
                this->communicate_shared_adof_ids(
                        tAdofListofTypes,
                        tDiscretizationIndexPerTypeAndTime,
                        tListSharedAdofIds,
                        tListSharedAdofPos );

                // Loop over all adofs determine the total number and the number of owned ones
                for ( uint Ik = 0; Ik < tNumTypeTimeLevels; Ik++ )
                {
                    // Set the Id of the shared adofs
                    for ( uint Ij = 0; Ij < tListSharedAdofIds( Ik ).numel(); Ij++ )
                    {
                        tAdofListofTypes( Ik )( tListSharedAdofPos( Ik )( Ij ) )->set_adof_id( tListSharedAdofIds( Ik )( Ij ) );
                    }
                }
            }

            // Tell pdofs to get adof Ids
            for ( uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
            {
                // all pdofs of this pdof host will ask for their adof Ids
                mPdofHostList( Ij )->get_adofs_ids();

                // create unique adof list for this pdof host
                // mPdofHostList(Ij)->create_unique_adof_list();
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Dof_Manager::set_pdof_t_matrix()
        {
            // Get number of pdof hosts in pdof host list
            uint tNumPdofHosts = mPdofHostList.size();

            // Tell pdofs to get the corresponding T matrix
            for ( uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
            {
                // all pdofs of this pdof host will ask for their t matrix
                mPdofHostList( Ij )->set_t_matrix( mUseHMR, mModelSolverInterface );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        Matrix< DDSMat >
        Dof_Manager::get_local_adof_ids()
        {
            Matrix< DDSMat > tLocalAdofIds( mNumOwnedAdofs, 1, -1 );

            // check for empty matrix
            if ( mNumOwnedAdofs > 0 )
            {
                uint tCounter = 0;

                for ( uint Ij = 0; Ij < mAdofListOwned.size(); Ij++ )
                {
                    for ( uint Ik = 0; Ik < mAdofListOwned( Ij ).size(); Ik++ )
                    {
                        tLocalAdofIds( tCounter++ ) = mAdofListOwned( Ij )( Ik )->get_adof_id();
                    }
                }

                MORIS_ASSERT( tLocalAdofIds.min() != -1,
                        "Dof_Manager::get_local_adof_ids(): Adof Id list not initialized correctly " );
            }

            return tLocalAdofIds;
        }

        //-----------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        Dof_Manager::get_local_adof_ids( const Vector< enum Dof_Type >& aListOfDofTypes )
        {
            Vector< Vector< Adof* > > tOwnedAdofs = this->get_owned_adofs();

            uint tCounter = 0;

            for ( uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                sint tPdofIndex = this->get_pdof_index_for_type( aListOfDofTypes( Ik ) );

                uint tTimeLevelOffsetForDofType = mTimeLevelOffsets( tPdofIndex );

                uint tTimeLevelPerType = mPdofHostTimeLevelList( tPdofIndex );

                for ( uint Ij = 0; Ij < tTimeLevelPerType; Ij++ )
                {
                    tCounter += tOwnedAdofs( tTimeLevelOffsetForDofType + Ij ).size();
                }
            }

            Matrix< DDSMat > tLocalAdofIds( tCounter, 1, -1 );

            tCounter = 0;

            for ( uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                sint tPdofIndex = this->get_pdof_index_for_type( aListOfDofTypes( Ik ) );

                uint tTimeLevelOffsetForDofType = mTimeLevelOffsets( tPdofIndex );

                uint tTimeLevelPerType = mPdofHostTimeLevelList( tPdofIndex );

                for ( uint Ij = 0; Ij < tTimeLevelPerType; Ij++ )
                {
                    for ( uint Ia = 0; Ia < tOwnedAdofs( tTimeLevelOffsetForDofType + Ij ).size(); Ia++ )
                    {
                        tLocalAdofIds( tCounter++ ) = tOwnedAdofs( tTimeLevelOffsetForDofType + Ij )( Ia )->get_adof_id();
                    }
                }
            }

            Matrix< DDSMat > tLocalAdofIdsSorted;

            sort( tLocalAdofIds, tLocalAdofIdsSorted );

            return tLocalAdofIdsSorted;
        }

        //-----------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        Dof_Manager::get_local_overlapping_adof_ids()
        {
            Matrix< DDSMat > tLocalAdofIds( mAdofList.size(), 1, -1 );

            // check for empty matrix
            if ( mAdofList.size() > 0 )
            {
                for ( uint Ij = 0; Ij < mAdofList.size(); Ij++ )
                {
                    tLocalAdofIds( Ij, 0 ) = mAdofList( Ij )->get_adof_id();
                }

                MORIS_ASSERT( tLocalAdofIds.min() != -1,
                        "Dof_Manager::get_local_overlapping_adof_ids(): Overlapping Adof Id list not initialized correctly " );
            }
            return tLocalAdofIds;
        }

        //-----------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        Dof_Manager::get_local_overlapping_adof_ids( const Vector< enum MSI::Dof_Type >& aListOfDofTypes )
        {
            uint tCounter = 0;

            for ( uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                sint tPdofIndex = this->get_pdof_index_for_type( aListOfDofTypes( Ik ) );

                uint tTimeLevelOffsetForDofType = mTimeLevelOffsets( tPdofIndex );

                uint tTimeLevelPerType = mPdofHostTimeLevelList( tPdofIndex );

                for ( uint Ij = 0; Ij < tTimeLevelPerType; Ij++ )
                {
                    tCounter += mAdofListOwnedAndShared( tTimeLevelOffsetForDofType + Ij ).size();
                }
            }

            Matrix< DDSMat > tLocalAdofIds( tCounter, 1, -1 );

            tCounter = 0;

            for ( uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                sint tPdofIndex = this->get_pdof_index_for_type( aListOfDofTypes( Ik ) );

                uint tTimeLevelOffsetForDofType = mTimeLevelOffsets( tPdofIndex );

                uint tTimeLevelPerType = mPdofHostTimeLevelList( tPdofIndex );

                for ( uint Ij = 0; Ij < tTimeLevelPerType; Ij++ )
                {
                    for ( uint Ia = 0; Ia < mAdofListOwnedAndShared( tTimeLevelOffsetForDofType + Ij ).size(); Ia++ )
                    {
                        tLocalAdofIds( tCounter++ ) = mAdofListOwnedAndShared( tTimeLevelOffsetForDofType + Ij )( Ia )->get_adof_id();
                    }
                }
            }

            return tLocalAdofIds;
        }

        //-----------------------------------------------------------------------------------------------------------

        moris::Matrix< DDSMat >
        Dof_Manager::get_unique_adof_mesh_indices()
        {
            sint tMaxNumAdofMeshInd = mModelSolverInterface->get_max_adof_index();

            moris::Matrix< DDSMat > tAdofMeshIndexExists( tMaxNumAdofMeshInd, 1, -1 );
            moris::Matrix< DDSMat > tAdofMeshIndex( tMaxNumAdofMeshInd, 1, -1 );
            uint                    tIndexCounter = 0;

            // Loop over all dof types. Determine the adof mesh index of these dofs
            for ( uint Ij = 0; Ij < mPdofTypeList.size(); Ij++ )
            {
                // Ask for adof order for this dof type
                uint tAdofMeshInd = (uint)mModelSolverInterface->get_adof_index_for_type( Ij );

                if ( tAdofMeshIndexExists( tAdofMeshInd ) == -1 )
                {
                    tAdofMeshIndexExists( tAdofMeshInd ) = 1;

                    tAdofMeshIndex( tIndexCounter++ ) = tAdofMeshInd;
                }
            }
            tAdofMeshIndex.resize( tIndexCounter, 1 );

            moris::Matrix< DDSMat > tAdofMeshIndexSorted;

            sort( tAdofMeshIndex, tAdofMeshIndexSorted );

            return tAdofMeshIndexSorted;
        }

        //-----------------------------------------------------------------------------------------------------------
        // this function is for HMR use only. It creates a map between MSI adof inds and HMR adof inds

        Matrix< DDUMat >
        Dof_Manager::get_adof_ind_map()
        {
            // Get length of adof list
            uint tAdofListSize = mAdofList.size();

            // Set size of adof list
            Matrix< DDUMat > tAdofIndMap( tAdofListSize, 1 );

            // Set external is of adof in mat
            for ( uint Ik = 0; Ik < tAdofListSize; Ik++ )
            {
                tAdofIndMap( Ik ) = mAdofList( Ik )->get_adof_external_ind();
            }

            return tAdofIndMap;
        }
    }    // namespace MSI
}    // namespace moris
