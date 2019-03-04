/*
 * cl_Dof_Manager.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_FEM_Node_Base.hpp"
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
        for ( moris::uint Ik = 0; Ik < mPdofHostList.size(); Ik++ )
        {
                delete mPdofHostList( Ik );
        }
        for ( moris::uint Ik = 0; Ik < mAdofList.size(); Ik++ )
        {
                delete mAdofList( Ik );
        }
    }

//-----------------------------------------------------------------------------------------------------------
    moris::uint Dof_Manager::initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // Ask how many equation objects
        moris::uint tNumEqnObj = aListEqnObj.size();
        moris::uint tMaxNumPdofHosts = 0;

        //loop over all equation objects, asking for their number of pdof hosts
        for ( moris::uint Ii=0; Ii < tNumEqnObj; Ii++ )
        {
            tMaxNumPdofHosts = std::max( tMaxNumPdofHosts, aListEqnObj( Ii )->get_max_pdof_hosts_ind() );
        }

        MORIS_ERROR( tMaxNumPdofHosts < MORIS_SINT_MAX, "Invalid Index in PDOF host" );

        // Add 1 because c++ is 0 based
        tMaxNumPdofHosts = tMaxNumPdofHosts + 1;

        return tMaxNumPdofHosts;
    }

//-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::initialize_pdof_type_list( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // Reserve of temporary pdof type list
        moris::Cell< enum Dof_Type > tTemporaryPdofTypeList;
        tTemporaryPdofTypeList.reserve( static_cast< int >( Dof_Type::END_ENUM ) + 1 );

        Matrix< DDUMat > tListToCheckIfEnumExist( (static_cast< int >(Dof_Type::END_ENUM) + 1), 1, 0 );

        // Get number of equation objects
        moris::uint tNumEquationObjects = aListEqnObj.size();

        //loop over all equation objects, asking for their pdof types
        for ( moris::uint Ii=0; Ii < tNumEquationObjects; Ii++ )
        {
            // Create temporary dof type list
            moris::Cell< enum Dof_Type > tDofType;

            // Ask equation object for its dof types
            aListEqnObj( Ii )->get_dof_types( tDofType );

            // Loop over all dof types
            for ( moris::uint Ik=0; Ik < tDofType.size(); Ik++ )
            {
                // Set 1 at position of the enum value
                tListToCheckIfEnumExist( static_cast< int >(tDofType( Ik )) ,0 ) = 1;
            }
        }

        // Loop over tListToCheckIfEnumExist. Add Ddof type to list if corresponding value = 1
        for ( moris::uint Ij=0; Ij < tListToCheckIfEnumExist.length(); Ij++ )
        {
            if ( tListToCheckIfEnumExist(Ij , 0) == 1)
            {
                tTemporaryPdofTypeList.push_back( ( Dof_Type ) Ij );
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
    void Dof_Manager::communicate_dof_types( moris::Cell< enum Dof_Type > & aPdofTypeList )
    {
        // Get processor size
        int tSize = par_size();

        // Get number of local dof types
        moris::sint tNumLocalDofTypes = aPdofTypeList.size();

        // Variable for maximal possible global dof types
        moris::sint tNumMaxGlobalDofTypes;

        // Get number of global dof types
        sum_all( tNumLocalDofTypes, tNumMaxGlobalDofTypes );

        if ( par_rank() == 0 )
        {
            // Set size of of pdof type list = number of global types
            mPdofTypeList.resize( tNumMaxGlobalDofTypes );
        }

        // Create list containing the number of local dof types
        moris::Cell < moris::sint > tNumLocalDofTypesList ( tSize );

        // Insert number of local dof types into list containing the number of local dof types
        MPI_Allgather( &tNumLocalDofTypes, 1, MPI_UNSIGNED, (tNumLocalDofTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD );

        // Create list containing the offsets of the local dof types in relation to processor 0
        moris::Cell< moris::sint > tDofTypeOffset( tSize, 0 );

        // Fill the list with the corresponding offsets
        for ( int Ip = 1; Ip < tSize; ++Ip )
        {
            tDofTypeOffset( Ip ) = tDofTypeOffset( Ip-1 ) + tNumLocalDofTypesList( Ip-1 );
        }

        // Assemble list containing all used dof types. Dof types are not unique
        MPI_Gatherv( ((aPdofTypeList.data()).data()),
                       tNumLocalDofTypes,
                       MPI_UNSIGNED,
                       (mPdofTypeList.data()).data(),
                       (tNumLocalDofTypesList.data()).data(),
                       (tDofTypeOffset.data()).data(),
                       MPI_UNSIGNED,
                       0,
                       MPI_COMM_WORLD );

        // Temporary variable for mPdofTypeList size
        moris::uint tPdofTypeListSize;

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
        MPI_Bcast( & tPdofTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // Resize mPdofTypeList on all processors
        mPdofTypeList.resize( tPdofTypeListSize );

        // Bcast unique mPdofTypeList to all processors
        MPI_Bcast( (mPdofTypeList.data()).data(), mPdofTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::create_dof_type_map()
    {
        //Get number of unique adofs of this equation object
        moris::uint tNumUniquePdofTypes = mPdofTypeList.size();

        // Get maximal dof type enum number
        moris::sint tMaxDofTypeEnumNumber = 0;

        // Loop over all pdof types to get the highest enum index
        for ( moris::uint Ii = 0; Ii < tNumUniquePdofTypes; Ii++ )
        {
            tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mPdofTypeList( Ii ) ) );
        }

        // +1 because c++ is 0 based
        tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;

        // Set size of maping matrix
        mPdofTypeMap       .set_size( tMaxDofTypeEnumNumber, 1, -1 );

        // Loop over all pdof types to create the mapping matrix
        for ( moris::uint Ii = 0; Ii < tNumUniquePdofTypes; Ii++ )
        {
            mPdofTypeMap( static_cast< int >( mPdofTypeList( Ii ) ), 0 ) = Ii;
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj )
     {
         moris::uint tNumUsedDofTypes = mPdofTypeList.size();
         moris::uint tNumEquationObj  = aListEqnObj.size();
         moris::uint tMaxNumPdofHosts = this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );

         // reserve the maximal needed amount of memory for a temporary PdofHostList
         moris::Cell < Pdof_Host* > tPdofHostList( tMaxNumPdofHosts, nullptr );

         // Loop over all equation objects. Ask them to create their pdof hosts. Pdof hosts are stored in tPdofHostList
         for ( moris::uint Ii=0; Ii < tNumEquationObj; Ii++ )
         {
             aListEqnObj( Ii )->create_my_pdof_hosts( tNumUsedDofTypes, mPdofTypeMap, tPdofHostList );
         }

         // Determine number of Pdof Hosts
         moris::uint tNumPdofHosts = 0;
         for ( moris::uint Ik = 0; Ik < tMaxNumPdofHosts; Ik++ )
         {
             // If pointer in temporary pdof host list exists add one to number of pdof hosts
             if ( tPdofHostList( Ik ) != NULL)
             {
                  tNumPdofHosts = tNumPdofHosts + 1;
             }
         }

         // Set size of List containing all pdof hosts
         mPdofHostList.resize( tNumPdofHosts );

         moris::uint counter = 0;
         // add pointers to pdof hosts into list of pdof hosts
         for ( moris::uint Ij = 0; Ij < tMaxNumPdofHosts; Ij++ )
         {
             // If pointer in temporary pdof host list exists add one to number of pdof hosts
             if ( tPdofHostList( Ij ) != NULL)
             {
                 mPdofHostList( counter ) = tPdofHostList( Ij );
                 counter = counter + 1;
             }
         }
     }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::initialize_pdof_host_time_level_list()
    {
        moris::uint tNumPdofTypes = mPdofTypeList.size();
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Set size of the list of pdof host time levels
        Matrix< DDUMat > tPdofHostTimeLevelList( tNumPdofTypes, 1, 0 );

        // Loop over pdof types
        for ( moris::uint Ik = 0; Ik < tNumPdofTypes; Ik++ )
        {
            // loop over pdof hosts
            for ( moris::uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
            {
                moris::uint tNumTimeLevels = mPdofHostList( Ii )->get_num_time_levels_of_type( Ik );

                if ( tNumTimeLevels != 0 )
                {
                    // Get the number of time levels per pdof type
                    tPdofHostTimeLevelList( Ik, 0 ) = tNumTimeLevels;

                    // Break loop if number of time levels for this pdof type is set.
                    break;
                }
            }
        }

        this->communicate_time_list( tPdofHostTimeLevelList );
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::communicate_time_list( Matrix< DDUMat > & aTimeLevelList )
    {
        mPdofHostTimeLevelList.set_size( mPdofTypeList.size(), 1, 0 );

        // Get the maximal value for all time list values
        for ( moris::uint Ii = 0; Ii < aTimeLevelList.length(); Ii++ )
        {
            max_all( aTimeLevelList( Ii, 0 ), mPdofHostTimeLevelList( Ii, 0 ) );
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::communicate_check_if_owned_adof_exists( moris::Cell< moris::Cell < Adof * > > & aAdofListofTypes )
    {
        // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
        Matrix< DDSMat > tCommTableMap ( mCommTable.max() + 1, 1, -1);

        moris::uint tNumCommProcs = mCommTable.length();

        // Loop over communication table to fill the communication table map
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            tCommTableMap( mCommTable( Ik, 0 ), 0 ) = Ik;
        }

        // Loop over all different adof types and time in this temporary list
        for ( moris::uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
        {
            moris::Cell< Matrix< DDUMat > > tSharedAdofPosGlobal( tNumCommProcs );

            // Set Mat to store number of shared adofs per processor
            Matrix< DDUMat > tNumSharedAdofsPerProc( tNumCommProcs, 1, 0 );

            // Loop over adofs per type/time. Count number of adofs per proc which have to be communicated
            for ( moris::uint Ib = 0; Ib < aAdofListofTypes( Ij ).size(); Ib++ )
            {
                // Check if adof at this position is not NULL
                if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    // Check if owning processor is this processor
                    if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                    {
                        // get owning procssor
                        moris::moris_id tProcID = aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                        moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                        MORIS_ASSERT( tProcIdPos != -1, "Dof_Manager::communicate_check_if_owned_adof_exists: Map returns proc rank -1. Check communication table");

                        // Add +1 to the processor number of shared dofs per processor
                        tNumSharedAdofsPerProc( tProcIdPos, 0) = tNumSharedAdofsPerProc( tProcIdPos, 0) + 1;
                    }
                }
            }

            // Set size of the moris::Mats in the Cell
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                if ( tNumSharedAdofsPerProc( Ik, 0 ) != 0 )
                {
                    tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                }
            }

            // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
            Matrix< DDUMat > tShredAdofPosPerProc( tNumCommProcs, 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ia = 0; Ia < aAdofListofTypes( Ij ).size(); Ia++ )
            {
                // Check if adof at this position is not NULL
                if ( aAdofListofTypes( Ij )( Ia ) != NULL )
                {
                    // Check if owning processor is this processor
                    if ( aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                    {
                        // Get owning procssor
                        moris::uint tProcID = aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor();

                        moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                        // Add owning procesor id to moris::Mat
                        tSharedAdofPosGlobal( tProcIdPos )( tShredAdofPosPerProc( tProcIdPos, 0 ), 0 ) = aAdofListofTypes( Ij )( Ia )->get_adof_external_id();

                        tShredAdofPosPerProc( tProcIdPos, 0 ) = tShredAdofPosPerProc( tProcIdPos, 0 ) + 1;
                    }
                }
            }

            // receiving list
            moris::Cell< Matrix< DDUMat > > tMatsToReceive;

            barrier();

            // Communicate position of shared adofs to the owning processor
            communicate_mats( mCommTable,
                              tSharedAdofPosGlobal,
                              tMatsToReceive );

            // Loop over all Mats set dummy owned adofs
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).length(); Ii++ )
                {
                    // Get owned adof Id
                    moris::uint tLocalAdofInd = mAdofGlobaltoLocalMap->find( tMatsToReceive( Ik )( Ii ) );

                    if ( aAdofListofTypes( Ij )( tLocalAdofInd ) == NULL )
                    {
                        //std::cout << "Invalid DOF Ownership. : " << par_rank() << " " << tLocalAdofInd << std::endl;
                        aAdofListofTypes( Ij )( tLocalAdofInd ) = new Adof();
                        aAdofListofTypes( Ij )( tLocalAdofInd )->set_adof_owning_processor( par_rank() );

                        // Set external adof ind. Used for HMR ordering
                        aAdofListofTypes( Ij )( tLocalAdofInd )->set_adof_external_ind( tLocalAdofInd );
                    }

//                    if ( aAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) == NULL )
//                    {
//                        aAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) = new Adof();
//                        aAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) )->set_adof_owning_processor( par_rank() );
//                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    moris::uint Dof_Manager::communicate_adof_offsets( const moris::uint & aNumOwnedAdofs )
    {
        // Get list containing the number of owned adofs of each processor
        Matrix< DDUMat > tNumOwnedAdofsList;
        comm_gather_and_broadcast( aNumOwnedAdofs, tNumOwnedAdofsList );

        Matrix< DDUMat > tOwnedAdofsOffsetList( tNumOwnedAdofsList.length(), 1, 0 );

        // Loop over all entries to create the offsets. Starting with 1
        for ( moris::uint Ij = 1; Ij < tOwnedAdofsOffsetList.length(); Ij++ )
        {
            // Add the number of owned adofs of the previous processor to the offset of the previous processor
            tOwnedAdofsOffsetList( Ij, 0 ) = tOwnedAdofsOffsetList( Ij-1, 0 ) + tNumOwnedAdofsList( Ij-1, 0 );
        }

        return tOwnedAdofsOffsetList( par_rank(), 0);
    }

    //----------------------------------------------------------------------------------------------------------
    void Dof_Manager::communicate_shared_adof_ids(const moris::Cell< moris::Cell < Adof * > > & aAdofListofTypes,
                                                        Matrix< DDUMat >             & aListSharedAdofIds,
                                                        Matrix< DDUMat >             & aListSharedAdofPos)
    {
        // Build communication table map to determine the right position for each processor rank. +1 because c++ is 0 based
        Matrix< DDSMat > tCommTableMap ( mCommTable.max() + 1, 1, -1);

        moris::uint tNumCommProcs = mCommTable.length();

        // Loop over communication table to fill the communication table map
        for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
        {
            tCommTableMap( mCommTable( Ik, 0 ), 0 ) = Ik;
        }

        moris::uint tCounter = 0;
        moris::uint tAdofPosCounter = 0;

        // Loop over all different adof types and time in this temporary list
        for ( moris::uint Ij = 0; Ij < aAdofListofTypes.size(); Ij++ )
        {
            moris::Cell< Matrix< DDUMat > > tSharedAdofPosGlobal( tNumCommProcs );
            moris::Cell< Matrix< DDUMat > > tSharedAdofPosLocal( tNumCommProcs );

            // Set Mat to store number of shared adofs per processor
            Matrix< DDUMat > tNumSharedAdofsPerProc( tNumCommProcs, 1, 0 );

            // Loop over adofs per type/time
            for ( moris::uint Ib = 0; Ib < aAdofListofTypes( Ij ).size(); Ib++ )
            {
                // Check if adof at this position is not NULL
                if ( aAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    // Check if owning processor is this processor
                    if ( aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                    {
                        // get owning procssor
                        moris::moris_id tProcID = aAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                        moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                        // Add +1 to the processor number of shared dofs per processor
                        tNumSharedAdofsPerProc( tProcIdPos, 0) = tNumSharedAdofsPerProc( tProcIdPos, 0) + 1;
                    }
                }
            }

            // Set size of the moris::Mats in the Cell
            for ( moris::uint Ik = 0; Ik < tNumCommProcs; Ik++ )
            {
                if ( tNumSharedAdofsPerProc( Ik, 0 ) != 0 )
                {
                    tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                    tSharedAdofPosLocal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                }
            }

            // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
            Matrix< DDUMat > tShredAdofPosPerProc( tNumCommProcs, 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ia = 0; Ia < aAdofListofTypes( Ij ).size(); Ia++ )
            {
                // Check if adof at this position is not NULL
                if ( aAdofListofTypes( Ij )( Ia ) != NULL )
                {
                    // Check if owning processor is this processor
                    if ( aAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                    {
                        // Get owning procssor
                        moris::moris_id tProcID = mAdofList( tCounter )->get_adof_owning_processor();

                        moris::sint tProcIdPos = tCommTableMap( tProcID, 0 );

                        // Add owning procesor id to moris::Mat
                        tSharedAdofPosGlobal( tProcIdPos )( tShredAdofPosPerProc( tProcIdPos, 0 ), 0 ) = mAdofList( tCounter )->get_adof_external_id();

                        // Add adof position to Mat
                        tSharedAdofPosLocal( tProcIdPos ) ( tShredAdofPosPerProc( tProcIdPos, 0 ), 0 ) = tCounter;

                        tShredAdofPosPerProc( tProcIdPos, 0 ) = tShredAdofPosPerProc( tProcIdPos, 0 ) + 1;
                    }
                    tCounter = tCounter + 1;
                }
            }

            // receiving list
            moris::Cell< Matrix< DDUMat > > tMatsToReceive;

            barrier();

            // Communicate position of shared adofs to the owning processor
            communicate_mats( mCommTable,
                              tSharedAdofPosGlobal,
                              tMatsToReceive );

            // Create List of Mats containing the shared node Ids
            moris::Cell< Matrix< DDUMat > > tSharesAdofIdList( tNumCommProcs );

            // Loop over all Mats setting the size
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                tSharesAdofIdList( Ik ).set_size( tMatsToReceive( Ik ).length(), 1);
            }

            // Loop over all received positions and get the adof id of the owning adof
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).length(); Ii++ )
                {
                    // Get owned adof Id
                    moris::uint tLocalAdofInd = mAdofGlobaltoLocalMap->find( tMatsToReceive( Ik )( Ii ) );

                    MORIS_ASSERT( ( aAdofListofTypes( Ij )( tLocalAdofInd )->get_adof_owning_processor() ) == par_rank(), "Dof_Manager::communicate_shared_adof_ids: Adof not owned by this processor");


                    tSharesAdofIdList( Ik )( Ii, 0 ) = ( aAdofListofTypes( Ij )( tLocalAdofInd ) )->get_adof_id();

                    //tSharesAdofIdList( Ik )( Ii, 0 ) = ( aAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) )->get_adof_id();
                }
            }

            moris::Cell< Matrix< DDUMat > > tMatsToReceive2;

            barrier();

            // Communicate owned adof Id back to the processor with the shared adof
            communicate_mats( mCommTable,
                              tSharesAdofIdList,
                              tMatsToReceive2 );

            // assemble Ids in list of shared adof ids and assemble the corresponding postions
            for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
            {
                if( tMatsToReceive2( Ik ).length() >= 1)
                {
                    aListSharedAdofIds ( {tAdofPosCounter, tAdofPosCounter + tMatsToReceive2( Ik ).length() -1 }, { 0, 0 } ) = tMatsToReceive2( Ik ).matrix_data();
                    aListSharedAdofPos ( {tAdofPosCounter, tAdofPosCounter +  tSharedAdofPosLocal( Ik ).length() -1 }, { 0, 0 } ) = tSharedAdofPosLocal( Ik ).matrix_data();

                    tAdofPosCounter =tAdofPosCounter + tMatsToReceive2( Ik ).length();
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::get_max_adof_ind( moris::sint & aMaxAdofInd )
    {
        // Get max entry of node adof if pdof host list exists
        if ( mNumMaxAdofs == -1 )
        {
            moris::uint tNumPdofHosts  = mPdofHostList.size();

            moris::Matrix< DDSMat> tAdofOrderExists( 3, 1, -1 );
            moris::Matrix< DDSMat> tAdofOrders( 3, 1, -1 );
            moris::uint tOrderCounter = 0;

            // Loop over all dof types. Determine the de adof orders of these dof types and put them into a list. can be 1,2,3
            for ( moris::uint Ij = 0; Ij < mPdofTypeList.size() ; Ij++ )
            {
                // Ask for adof order for this dof type
                moris::uint tAdofOrder = ( moris::uint ) mModelSolverInterface->get_adof_order_for_type( Ij );

                if( tAdofOrderExists( tAdofOrder-1, 0 ) == -1)
                {
                     tAdofOrderExists( tAdofOrder-1, 0 ) = 1;

                     tAdofOrders( tOrderCounter++ ) = tAdofOrder;
                }
            }
            tAdofOrders.resize( tOrderCounter, 1 );

            if ( tNumPdofHosts != 0 )
            {
                for ( moris::uint Ia = 0; Ia < tAdofOrders.length() ; Ia++ )
                {
                    moris::uint tAdofOrder = tAdofOrders( Ia, 0 );

                    for ( moris::uint Ik = 0; Ik < tNumPdofHosts; Ik++ )
                    {
                        moris::fem::Node_Base * tNode = mPdofHostList( Ik )->get_node_obj_ptr();

                        aMaxAdofInd = std::max( aMaxAdofInd, ( tNode->get_adof_indices( tAdofOrder ) ).max() );
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
        //        else if ( mUseHMR == false )
        //        {
        //            aMaxAdofInd = mPdofHostList.size();
        //        }
        else { MORIS_ERROR( false, "MSI::Dof_Manager: Check number of adofs"); }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::create_adofs()
    {
        this->initialize_pdof_host_time_level_list();

        Matrix< DDUMat > tTimeLevelOffsets( mPdofHostTimeLevelList.length(), 1, 0);
        for ( moris::uint Ik = 1; Ik < mPdofHostTimeLevelList.length(); Ik++ )
        {
            tTimeLevelOffsets( Ik, 0 ) = tTimeLevelOffsets( Ik-1, 0 ) + mPdofHostTimeLevelList( Ik-1, 0 );
        }

        // Get number of pdoftypes and size of pdof host list
        moris::uint tNumPdofHosts  = mPdofHostList.size();
        moris::uint tNumTimeLevels = sum( mPdofHostTimeLevelList );
        moris::sint tMaxNodeAdofId = -1;

        this->get_max_adof_ind( tMaxNodeAdofId );

        // Create temporary moris::Cell containing lists of temporary adofs
        moris::Cell<moris::Cell < Adof * > > tAdofListofTypes( tNumTimeLevels );

        for ( moris::uint Ik = 0; Ik < tNumTimeLevels; Ik++ )
        {
            tAdofListofTypes( Ik ).resize( tMaxNodeAdofId, nullptr );
        }

        // Loop over all pdof hosts and get the adofs
        for ( moris::uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
        {
            mPdofHostList( Ii )->get_adofs( tTimeLevelOffsets, tAdofListofTypes, mModelSolverInterface, mUseHMR );
        }

        // Check if shared adof exists
        if ( !(par_size() <= 1) )
        {
            this->communicate_check_if_owned_adof_exists( tAdofListofTypes );
        }

        // Determine number of owned and shared adofs
        moris::uint tNumAdofs = 0;
        moris::uint tNumOwnedAdofs = 0;
        moris::uint tNumSharedAdofs = 0;

        // Multigrid type and time identifier;
        moris::sint tAdofTypeTimeIdentifier = 0;
        mTypeTimeIndentifierToTypeMap.set_size( tAdofListofTypes.size(), 1 , -1 );

        // Loop over all adofs determine the total number and the number of owned ones
        for ( moris::uint Ik = 0; Ik < tAdofListofTypes.size(); Ik++ )
        {
            for ( moris::uint Ia = 0; Ia < tAdofListofTypes( Ik ).size(); Ia++ )
            {
                // If pointer in temporary adof list exists. Add one to number of owned adofs
                if ( ( tAdofListofTypes( Ik )( Ia ) != NULL ) )
                {
                    tNumAdofs = tNumAdofs + 1;

                    // If owning processor equals this processor then its an owned adof
                    if ( ( tAdofListofTypes( Ik )( Ia )->get_adof_owning_processor() == par_rank() ) )
                    {
                        tNumOwnedAdofs = tNumOwnedAdofs + 1;
                    }

                    // Add type/time identifier to Adof
                    tAdofListofTypes( Ik )( Ia )->set_adof_type_time_identifier( tAdofTypeTimeIdentifier );
                }
            }

            // Multigrid Type time identifier to type map
            moris::uint tCounterTypeTime = 1;
            if ( tCounterTypeTime < tTimeLevelOffsets.length() )
            {
                if ( Ik < tTimeLevelOffsets( tCounterTypeTime, 0 ))
                {
                    mTypeTimeIndentifierToTypeMap( Ik, 0 ) = tCounterTypeTime-1;
                }
                else if ( Ik == tTimeLevelOffsets( tCounterTypeTime, 0 ) )
                {
                    tCounterTypeTime++;
                    mTypeTimeIndentifierToTypeMap( Ik, 0 ) = tCounterTypeTime-1;
                }
                else
                {
                    MORIS_ASSERT( false, "Dof_Manager::create_adofs(), " );
                }
            }
            else
            {
                mTypeTimeIndentifierToTypeMap( Ik, 0 ) = tCounterTypeTime-1;
            }

            tAdofTypeTimeIdentifier++;
        }

        // Calculate number of shared adofs
        tNumSharedAdofs = tNumAdofs - tNumOwnedAdofs;

        // Get adof offset for this processor
        moris::uint tAdofOffset = this->communicate_adof_offsets( tNumOwnedAdofs );

        // Set size of List containing all adofs
        mAdofList.resize( tNumAdofs );
        mAdofListOwned.resize( tNumOwnedAdofs );
        //mAdofListShared.resize( tNumSharedAdofs );

        moris::uint tCounterOwned = tAdofOffset;
        //moris::uint tCounterShared = 0;
        moris::uint tCounter = 0;

        // loop over temporary adof list. Add pointers to adofs into list of adofs
        for ( moris::uint Ij = 0; Ij < tAdofListofTypes.size(); Ij++ )
        {
            for ( moris::uint Ib = 0; Ib < tAdofListofTypes( Ij ).size(); Ib++ )
            {
                // If pointer in temporary adofs list exists, add adof
                if ( tAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    if ( tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() == par_rank() )
                    {
                        // Add adof to owned adof list
                        mAdofListOwned( tCounterOwned - tAdofOffset ) = tAdofListofTypes( Ij )( Ib ) ;

                        // Set adof Id
                        mAdofListOwned( tCounterOwned - tAdofOffset )->set_adof_id( tCounterOwned );

                        tCounterOwned = tCounterOwned + 1;
                    }
//                    else
//                    {
//                        mAdofListShared( tCounterShared ) = tAdofListofTypes( Ij )( Ib );
//                        tCounterShared = tCounterShared + 1;
//                    }
                    // Add adof to adof list
                    mAdofList( tCounter ) = tAdofListofTypes( Ij )( Ib ) ;

                    tCounter = tCounter + 1;
                }
            }
        }

        // Create vector storing information about adof ids and position for communication
        Matrix< DDUMat > tListSharedAdofIds;
        Matrix< DDUMat > tListSharedAdofPos;

        if ( !(par_size() <= 1) )
        {
            tListSharedAdofIds.set_size( tNumSharedAdofs, 1 );
            tListSharedAdofPos.set_size( tNumSharedAdofs, 1 );

            // Communicate shared-owned adof ids
            this->communicate_shared_adof_ids( tAdofListofTypes, tListSharedAdofIds, tListSharedAdofPos );

            // Set the Id of the shared adofs
            for ( moris::uint Ij = 0; Ij < tListSharedAdofIds.length(); Ij++ )
            {
                mAdofList( tListSharedAdofPos( Ij ) )->set_adof_id( tListSharedAdofIds( Ij ) );
            }
        }

        // Tell pdofs to get adof Ids
        for ( moris::uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
        {
            // all pdofs of this pdof host will ask for their adof Ids
            mPdofHostList(Ij)->get_adofs_ids();

            // create unique adof list for this pdof host
            //mPdofHostList(Ij)->create_unique_adof_list();
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Dof_Manager::set_pdof_t_matrix()
    {
        // Get number of pdof hosts in pdof host list
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Tell pdofs to get t their t matrix
        for ( moris::uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
        {
            // all pdofs of this pdof host will ask for their t matrix
            mPdofHostList(Ij)->set_t_matrix( mUseHMR, mModelSolverInterface );
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    Matrix< DDSMat > Dof_Manager::get_local_adof_ids()
    {
        Matrix< DDSMat > tLocalAdofIds ( mAdofListOwned.size(), 1 );

        for ( moris::uint Ij = 0; Ij < mAdofListOwned.size(); Ij++ )
        {
            tLocalAdofIds( Ij, 0 ) = mAdofListOwned( Ij )->get_adof_id();
        }

        MORIS_ASSERT( tLocalAdofIds.min() != -1, "Dof_Manager::get_local_adof_ids(): Adof Id list not initialized correctly ");

        return tLocalAdofIds;
    }

    //-----------------------------------------------------------------------------------------------------------

    Matrix< DDSMat > Dof_Manager::get_local_adof_ids( const moris::Cell< enum Dof_Type > & aListOfDofTypes )
    {
        // Initialize counter
        moris::uint tCounterAdofIds = 0;//FIXME check if owned

        // Loop over all pdof hosts
        for ( moris::uint Ij = 0; Ij < mPdofHostList.size(); Ij++ )//FIXME check if owned
        {
            // Loop over all dof types
            for ( moris::uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )//FIXME check if owned
            {
                // Get dof type index
                moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aListOfDofTypes( Ik ) ) );

                // get number of time levels on this dof type
                moris::uint tTimeLevels = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex ).size();
                for ( moris::uint Ii = 0; Ii < tTimeLevels; Ii++ )
                {
                    // Get vector with adof ids for this pdof
                    Matrix< DDSMat > tAdofIds = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex )( Ii )->mAdofIds;

                    tCounterAdofIds =+ tAdofIds.length();
                }
            }
        }

        // Initialize
        Matrix< DDSMat > tLocalAdofIds ( tCounterAdofIds, 1, -1 );

        // Re-initialize counter
        tCounterAdofIds = 0;

        // Loop over all pdof hosts
        for ( moris::uint Ij = 0; Ij < mPdofHostList.size(); Ij++ )
        {
            // Loop over all dof types
            for ( moris::uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                // Get dof type index
                moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aListOfDofTypes( Ik ) ) );

                // get number of time levels on this dof type
                moris::uint tTimeLevels = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex ).size();
                for ( moris::uint Ii = 0; Ii < tTimeLevels; Ii++ )
                {
                    // Get vector with adof ids for this pdof
                    Matrix< DDSMat > tAdofIds = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex )( Ii )->mAdofIds;

                    // Add adof Ids to list
                    tLocalAdofIds( {tCounterAdofIds, tCounterAdofIds + tAdofIds.length() -1 }, { 0, 0 } ) = tAdofIds.matrix_data();

                    tCounterAdofIds =+ tAdofIds.length();
                }
            }
        }

        // make list unique
        Matrix< DDSMat > tLocalUniqueAdofIds;
        unique( tLocalAdofIds, tLocalUniqueAdofIds );

        MORIS_ASSERT( tLocalUniqueAdofIds.min() != -1, "Dof_Manager::get_local_adof_ids(): Adof Id list not initialized correctly ");

        return tLocalUniqueAdofIds;
    }

    //-----------------------------------------------------------------------------------------------------------
    Matrix< DDSMat > Dof_Manager::get_local_overlapping_adof_ids( const moris::Cell< enum Dof_Type > & aListOfDofTypes )
    {
        // Initialize counter
        moris::uint tCounterAdofIds = 0;

        // Loop over all pdof hosts
        for ( moris::uint Ij = 0; Ij < mPdofHostList.size(); Ij++ )
        {
            // Loop over all dof types
            for ( moris::uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                // Get dof type index
                moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aListOfDofTypes( Ik ) ) );

                // get number of time levels on this dof type
                moris::uint tTimeLevels = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex ).size();
                for ( moris::uint Ii = 0; Ii < tTimeLevels; Ii++ )
                {
                    // Get vector with adof ids for this pdof
                    Matrix< DDSMat > tAdofIds = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex )( Ii )->mAdofIds;

                    tCounterAdofIds =+ tAdofIds.length();
                }
            }
        }

        // Initialize
        Matrix< DDSMat > tLocalAdofIds ( tCounterAdofIds, 1, -1 );

        // Re-initialize counter
        tCounterAdofIds = 0;

        // Loop over all pdof hosts
        for ( moris::uint Ij = 0; Ij < mPdofHostList.size(); Ij++ )
        {
            // Loop over all dof types
            for ( moris::uint Ik = 0; Ik < aListOfDofTypes.size(); Ik++ )
            {
                // Get dof type index
                moris::sint tDofTypeIndex = mPdofTypeMap( static_cast< int >( aListOfDofTypes( Ik ) ) );

                // get number of time levels on this dof type
                moris::uint tTimeLevels = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex ).size();
                for ( moris::uint Ii = 0; Ii < tTimeLevels; Ii++ )
                {
                    // Get vector with adof ids for this pdof
                    Matrix< DDSMat > tAdofIds = mPdofHostList( Ik )->get_pdof_time_list( tDofTypeIndex )( Ii )->mAdofIds;

                    // Add adof Ids to list
                    tLocalAdofIds( {tCounterAdofIds, tCounterAdofIds + tAdofIds.length() -1 }, { 0, 0 } ) = tAdofIds.matrix_data();

                    tCounterAdofIds =+ tAdofIds.length();
                }
            }
        }

        // make list unique
        Matrix< DDSMat > tLocalUniqueAdofIds;
        unique( tLocalAdofIds, tLocalUniqueAdofIds );

        MORIS_ASSERT( tLocalUniqueAdofIds.min() != -1, "Dof_Manager::get_local_adof_ids(): Adof Id list not initialized correctly ");

        return tLocalUniqueAdofIds;
    }

    //-----------------------------------------------------------------------------------------------------------
    Matrix< DDSMat > Dof_Manager::get_local_overlapping_adof_ids()
    {
        Matrix< DDSMat > tLocalAdofIds ( mAdofList.size(), 1, -1 );

        for ( moris::uint Ij = 0; Ij < mAdofList.size(); Ij++ )
        {
            tLocalAdofIds( Ij, 0 ) = mAdofList( Ij )->get_adof_id();
        }

        MORIS_ASSERT( tLocalAdofIds.min() != -1, "Dof_Manager::get_local_overlapping_adof_ids(): Overlapping Adof Id list not initialized correctly ");

        return tLocalAdofIds;
    }

    //-----------------------------------------------------------------------------------------------------------
    moris::Matrix< DDSMat > Dof_Manager::get_unique_dof_type_orders()
    {
        moris::Matrix< DDSMat> tAdofOrderExists( 3, 1, -1 );
        moris::Matrix< DDSMat> tAdofOrders( 3, 1, -1 );
        moris::uint tOrderCounter = 0;

        // Loop over all dof types. Determine the de adof orders of these dof types and put them into a list. can be 1,2,3
        for ( moris::uint Ij = 0; Ij < mPdofTypeList.size() ; Ij++ )
        {
            // Ask for adof order for this dof type
            moris::uint tAdofOrder = ( moris::uint ) mModelSolverInterface->get_adof_order_for_type( Ij );

            if( tAdofOrderExists( tAdofOrder-1, 0 ) == -1)
            {
                 tAdofOrderExists( tAdofOrder-1, 0 ) = 1;

                 tAdofOrders( tOrderCounter++ ) = tAdofOrder;
            }
        }
        tAdofOrders.resize( tOrderCounter, 1 );

        moris::Matrix< DDSMat> tAdofOrdersSorted;

        sort( tAdofOrders, tAdofOrdersSorted );

        return tAdofOrdersSorted;
    }

    //-----------------------------------------------------------------------------------------------------------
    //this function is for HMR use only. It creates a map between MSI adof inds and HMR adof inds
    Matrix< DDUMat > Dof_Manager::get_adof_ind_map()
    {
        // Get length of adof list
        moris::uint tAdofListSize = mAdofList.size();

        // Set size of adof list
        Matrix< DDUMat > tAdofIndMap( tAdofListSize, 1 );

        // Set external is of adof in mat
        for ( moris::uint Ik = 0; Ik < tAdofListSize; Ik++ )
        {
            tAdofIndMap( Ik , 0 ) = mAdofList( Ik )->get_adof_external_ind();
        }

        return tAdofIndMap;
    }



}
}
