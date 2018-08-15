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
#include "cl_Equation_Object.hpp"

namespace moris
{
namespace MSI
{
class Pdof_Host;
class Dof_Manager
{
private:
    moris::Cell < Pdof_Host * > mPdofHostList;
    moris::Cell < Adof * > mAdofList;
    moris::Cell < Adof * > mAdofListOwned;

    moris::uint mMaxNumPdofHosts;

    // List containing all dof types. Perhaps make temporary in future.
    moris::Cell< enum Dof_Type > mPdofTypeList;

    moris::Mat< moris::sint > mPdofTypeMap;

public:
    Dof_Manager()
    {
        //this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );
    };

    Dof_Manager( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );

        this->initialize_pdof_type_list( aListEqnObj );

        this->initialize_pdof_host_list( aListEqnObj );

        this->create_adofs();

        this->set_pdof_t_matrix();

        for ( moris::uint Ii=0; Ii < aListEqnObj.size(); Ii++ )
        {
            aListEqnObj( Ii )->create_my_pdof_list();
            aListEqnObj( Ii )->create_my_list_of_adof_ids();

            aListEqnObj( Ii )->set_unique_adof_map();
        }
    };

    ~Dof_Manager()
    {};

//-----------------------------------------------------------------------------------------------------------
    void initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // Ask how many equation objects
        moris::uint tNumEqnObj = aListEqnObj.size();
        mMaxNumPdofHosts = 0;

        //loop over all equation objects, asking for their number of pdof hosts
        for ( moris::uint Ii=0; Ii < tNumEqnObj; Ii++ )
        {
            mMaxNumPdofHosts += aListEqnObj( Ii )->get_num_pdof_hosts();
        }
    };

    //-----------------------------------------------------------------------------------------------------------
    void initialize_pdof_type_list( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // reserve 512 slots  // perhaps reserve max number of enums
        moris::Cell< enum Dof_Type > tPdofTypeList2;
        tPdofTypeList2.reserve( 512 );                          //FIXME typdef

        //loop over all equation objects, asking for their pdof types
        for ( moris::uint Ii=0; Ii < aListEqnObj.size(); Ii++ )
        {
            // Create temporary dof type list
            moris::Cell< enum Dof_Type > tDofType;

            // Ask equation object for its dof types
            aListEqnObj( Ii )->get_dof_types( tDofType );

            //loop over all equation objects, asking for their pdof types
            for ( moris::uint Ik=0; Ik < tDofType.size(); Ik++ )
            {
                bool tDofTypeExists = false;

                // loop over all dof types of this equation object
                for ( moris::uint Ij=0; Ij < tPdofTypeList2.size(); Ij++ )
                {
                    // if dof type exists set tDofTypeExists to true
                    if( tPdofTypeList2( Ij ) == tDofType( Ik ) )
                    {
                        // Check if dof tupe exists twice
                        MORIS_ERROR( ! tDofTypeExists, "Pdof_Host::set_dof_type(): Dof type exists twice. This is not allowed");

                        tDofTypeExists = true;
                        // break;
                    }
                }

                // If dof type does not exist tDofTypeExists = false, then add dof type to dof type list
                if( ! tDofTypeExists )
                {
                    tPdofTypeList2.push_back( tDofType( Ik ) );
                }
            }
        }
        tPdofTypeList2.shrink_to_fit();

        this->communicate_dof_types( tPdofTypeList2 );

        this->create_dof_type_map();
    };

    //-----------------------------------------------------------------------------------------------------------
    void communicate_dof_types( moris::Cell< enum Dof_Type > & aPdofTypeList )
    {
        // Get number of local dof types
        moris::sint tNumLocalDofTypes = aPdofTypeList.size();

        // Variable for maximal possible global dof types
        moris::sint tNumMaxGlobalDofTypes;

        // Get number of global dof types
        MPI_Allreduce( &tNumLocalDofTypes, &tNumMaxGlobalDofTypes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

        // Get processor size and rank
        //int tRank = par_rank();
        int tSize = par_size();

        // Set size of of pdof type list = number of global types
        mPdofTypeList.resize( tNumMaxGlobalDofTypes );

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
        MPI_Allgatherv( ((aPdofTypeList.data()).data()),
                        tNumLocalDofTypes,
                        MPI_UNSIGNED,
                        (mPdofTypeList.data()).data(),
                        (tNumLocalDofTypesList.data()).data(),
                        (tDofTypeOffset.data()).data(),
                        MPI_UNSIGNED,
                        MPI_COMM_WORLD );

//        if(tRank == 0)
//        {
//            for ( moris::uint Ij=0; Ij < mPdofTypeList1.size(); Ij++ )
//            {
//            std::cout<<static_cast<int>(mPdofTypeList1(Ij))<<std::endl;
//            }
//        }

        // Sort this created list
        std::sort( (mPdofTypeList.data()).data(), (mPdofTypeList.data()).data() + mPdofTypeList.size() );

        // use std::unique and std::distance to create  list containing all used dof types. This list is unique
        auto last = std::unique( (mPdofTypeList.data()).data(), (mPdofTypeList.data()).data() + mPdofTypeList.size() );
        auto pos  = std::distance( (mPdofTypeList.data()).data(), last );

        mPdofTypeList.resize( pos );

//      if(tRank == 0)
//      {
//          for ( moris::uint Ij=0; Ij < mPdofTypeList.size(); Ij++ )
//          {
//          std::cout<<static_cast<int>(mPdofTypeList(Ij))<<std::endl;
//          }
//      }
    };

    //-----------------------------------------------------------------------------------------------------------
    void create_dof_type_map()
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
        mPdofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

        // Loop over all pdof types to create the mapping matrix
        for ( moris::uint Ii = 0; Ii < tNumUniquePdofTypes; Ii++ )
        {
            mPdofTypeMap( static_cast< int >( mPdofTypeList( Ii ) ), 0 ) = Ii;
        }
    };

//-----------------------------------------------------------------------------------------------------------
    void initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        moris::uint tNumUsedDofTypes = mPdofTypeList.size();
        moris::uint tNumEquationObj  = aListEqnObj.size();

        // reserve the maximal needed amount of memory for  a temporary PdofHostList
        moris::Cell < Pdof_Host* > tPdofHostList( mMaxNumPdofHosts, nullptr );

        // Loop over all equation objects. Ask them to create their pdof hosts. Pdof hosts are stored in tPdofHostList
        for ( moris::uint Ii=0; Ii < tNumEquationObj; Ii++ )
        {
            aListEqnObj( Ii )->create_my_pdof_hosts( tNumUsedDofTypes, mPdofTypeMap, tPdofHostList );
        }

        // Determine number of Pdof Hosts
        moris::uint tNumPdofHosts = 0;
        for ( moris::uint Ik = 0; Ik < mMaxNumPdofHosts; Ik++ )
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
        for ( moris::uint Ij = 0; Ij < mMaxNumPdofHosts; Ij++ )
        {
            // If pointer in temporary pdof host list exists add one to number of pdof hosts
            if ( tPdofHostList( Ij ) != NULL)
            {
                mPdofHostList( counter ) = tPdofHostList( Ij );
                counter = counter + 1;
            }
        }
    };

//-----------------------------------------------------------------------------------------------------------
    void create_adofs()
    {
        // Get number of pdoftypes
        moris::uint tNumPdofTypes = mPdofTypeList.size();

        // Get number of pdof hosts in pdof host list
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Get max entry of node adof
        moris::sint tMaxNodeAdofId = 0;
        for ( moris::uint Ik = 0; Ik < tNumPdofHosts; Ik++ )
        {
            tMaxNodeAdofId = std::max( tMaxNodeAdofId, mPdofHostList( Ik )->get_node_obj_ptr()->get_adof_ids().max() );
        }

        // Add one because c++ is 0 based. ==> List size has to be tMaxNodeAdofId + 1
        tMaxNodeAdofId = tMaxNodeAdofId +1;

        // Create temporary moris::Cell containing lists of temporary adofs
        moris::Cell<moris::Cell < Adof * > > tAdofListofTypes;
        tAdofListofTypes.resize( tNumPdofTypes );

        for ( moris::uint Ik = 0; Ik < tNumPdofTypes; Ik++ )
        {
            tAdofListofTypes( Ik ).resize( tMaxNodeAdofId, nullptr );
        }

        // Loop over all pdof hosts and get the adofs
        for ( moris::uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
        {
            mPdofHostList( Ii )->get_adofs( tAdofListofTypes );
        }

        if ( !(par_size() <= 1) )
        {
            this->communicate_check_if_owned_adof_exists( tAdofListofTypes );
        }

        // Determine number of owned and shared adofs
        moris::uint tNumAdofs = 0;
        moris::uint tNumOwnedAdofs = 0;
        moris::uint tNumSharedAdofs = 0;

        for ( moris::uint Ik = 0; Ik < tAdofListofTypes.size(); Ik++ )
        {
            for ( moris::uint Ia = 0; Ia < tAdofListofTypes( Ik ).size(); Ia++ )
            {
                // If pointer in temporary adof list exists. Add one to number of owned adofs
                if ( ( tAdofListofTypes( Ik )( Ia ) != NULL ) )
                {
                    tNumAdofs = tNumAdofs + 1;

                    if ( ( tAdofListofTypes( Ik )( Ia )->get_adof_owning_processor() == par_rank() ) )
                    {
                        tNumOwnedAdofs = tNumOwnedAdofs + 1;
                    }
                }
            }
        }

        // Calculate number of shared adofs
        tNumSharedAdofs = tNumAdofs - tNumOwnedAdofs;

        // Get adof offset for this processor
        moris::uint tAdofOffset = this->communicate_adof_offsets( tNumOwnedAdofs );

        // Set size of List containing all adofs
        mAdofList.resize( tNumAdofs );
        mAdofListOwned.resize( tNumOwnedAdofs );

        moris::uint tCounterOwned = tAdofOffset;
        moris::uint tCounter = 0; //FIXME

        // add pointers to adofs into list of adofs
        for ( moris::uint Ij = 0; Ij < tAdofListofTypes.size(); Ij++ )
        {
            for ( moris::uint Ib = 0; Ib < tAdofListofTypes( Ij ).size(); Ib++ )
            {
                // If pointer in temporary adofs list exists, add one to number of adofs
                if ( tAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    if (  tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() == par_rank() )
                    {
                        mAdofListOwned( tCounterOwned - tAdofOffset ) = tAdofListofTypes( Ij )( Ib ) ;

                        // Set adof Id. Fixme introduce offset for parallel
                        mAdofListOwned( tCounterOwned - tAdofOffset )->set_adof_id( tCounterOwned );

                        tCounterOwned = tCounterOwned + 1;
                    }
                    // FIXME ovverring pointer
                    mAdofList( tCounter ) = tAdofListofTypes( Ij )( Ib ) ;
                    //mAdofList( tCounter )->set_adof_id( tCounter );
                    tCounter = tCounter + 1;
                }
            }
        }

        moris::Mat< moris::uint > tListSharedAdofIds;
        moris::Mat< moris::uint > tListSharedAdofPos;

        if ( !(par_size() <= 1) )
        {
            tListSharedAdofIds.set_size( tNumSharedAdofs, 1 );
            tListSharedAdofPos.set_size( tNumSharedAdofs, 1 );
            this->communicate_shared_adof_ids( tAdofListofTypes, tListSharedAdofIds, tListSharedAdofPos );

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
            mPdofHostList(Ij)->create_unique_adof_list();
        }
    };

    //-----------------------------------------------------------------------------------------------------------
    moris::uint communicate_adof_offsets( const moris::uint & aNumOwnedAdofs )
    {
        // Get list containing the number of owned adofs of each processor
        moris::Mat< moris::uint > tNumOwnedAdofsList = comm_gather_and_broadcast( aNumOwnedAdofs );

        moris::Mat< moris::uint > tOwnedAdofsOffsetList( tNumOwnedAdofsList.length(), 1, 0 );

        // Loop over all entries to create the offsets. Starting with 1
        for ( moris::uint Ij = 1; Ij < tOwnedAdofsOffsetList.length(); Ij++ )
        {
            // Add the number of owned adofs of the previous processor to the offset of the previous processor
            tOwnedAdofsOffsetList( Ij, 0 ) = tOwnedAdofsOffsetList( Ij-1, 0 ) + tNumOwnedAdofsList( Ij-1, 0 );
        }

        return tOwnedAdofsOffsetList( par_rank(), 0);
    };


    //-----------------------------------------------------------------------------------------------------------
    void communicate_shared_adof_ids(const moris::Cell< moris::Cell < Adof * > > & tAdofListofTypes,
                                           moris::Mat< moris::uint > & aAaaa,
                                           moris::Mat< moris::uint > & aBbbb)
    {
        moris::Mat< moris::uint > aCommunicationList( par_size(), 1 );

        for ( moris::uint Ik = 0; Ik < aCommunicationList.length(); Ik++ )
        {
            aCommunicationList( Ik, 0 ) = Ik;
        }

        moris::uint tCounter = 0;
        moris::uint tAdofPosCounter = 0;

        // Loop over all different adof types in this temporary list
        for ( moris::uint Ij = 0; Ij < tAdofListofTypes.size(); Ij++ )
        {
            moris::Cell< moris::Mat< moris::uint > > tSharedAdofPosGlobal( par_size() );    //What to ask
            moris::Cell< moris::Mat< moris::uint > > tSharedAdofPosLocal( par_size() );     // where to place

            // Set Mat to store number of shared adofs per processor
            moris::Mat< moris::uint > tNumSharedAdofsPerProc( par_size(), 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ib = 0; Ib < tAdofListofTypes( Ij ).size(); Ib++ )
            {
                // Check if adof at this position is not NULL
                if ( tAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    // Check if owning processor is this processor
                    if (  tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                    {
                        // get owning procssor
                        moris::uint tProcID = tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                        // Add +1 to the processor number of shared dofs per processor
                        tNumSharedAdofsPerProc( tProcID, 0) = tNumSharedAdofsPerProc( tProcID, 0) + 1;
                    }
                }
            }

            // Set size of the moris::Mats in the Cell
            for ( moris::uint Ik = 0; Ik < par_size(); Ik++ )
            {
                if ( tNumSharedAdofsPerProc( Ik, 0 ) != 0 )
                {
                    tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                    tSharedAdofPosLocal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                }
            }

            // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
            moris::Mat< moris::uint > tShredAdofPosPerProc( par_size(), 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ia = 0; Ia < tAdofListofTypes( Ij ).size(); Ia++ )
            {
                // Check if adof at this position is not NULL
                if ( tAdofListofTypes( Ij )( Ia ) != NULL )
                {
                    // Check if owning processor is this processor
                    if (  tAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                    {
                        // Get owning procssor
                        moris::uint tProcId = mAdofList( tCounter )->get_adof_owning_processor();

                        // Add owning procesor id to moris::Mat
                        tSharedAdofPosGlobal( tProcId )( tShredAdofPosPerProc( tProcId, 0 ), 0 ) = mAdofList( tCounter )->get_adof_external_id();

                        // Add adof position to Mat
                        tSharedAdofPosLocal( tProcId ) ( tShredAdofPosPerProc( tProcId, 0 ), 0 ) = tCounter;

                        tShredAdofPosPerProc( tProcId, 0 ) = tShredAdofPosPerProc( tProcId, 0 ) + 1;
                    }
                    tCounter = tCounter + 1;
                }
            }

            // receiving list
            moris::Cell< moris::Mat< moris::uint > > tMatsToReceive;

            barrier();

            communicate_mats( aCommunicationList,
                              tSharedAdofPosGlobal,
                              tMatsToReceive );

            // Create List of Mats containing the shared node Ids
            moris::Cell< moris::Mat< moris::uint > > tSharesAdofIdList( par_size() );

            // Loop over all Mats setting the size
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                tSharesAdofIdList( Ik ).set_size( tMatsToReceive( Ik ).length(), 1);
            }

            // Loop iver all Mats
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).length(); Ii++ )
                {
                    // Get owned adof Id
                    tSharesAdofIdList( Ik )( Ii, 0 ) = (tAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) )->get_adof_id();
                }
            }

            moris::Cell< moris::Mat< moris::uint > > tMatsToReceive2;

            communicate_mats( aCommunicationList,
                              tSharesAdofIdList,
                              tMatsToReceive2 );

            // assemble in Mat
            for ( moris::uint Ik = 0; Ik < tMatsToReceive2.size(); Ik++ )
            {
                if( tMatsToReceive2( Ik ).length() >= 1)
                {
                    aAaaa ( {tAdofPosCounter, tAdofPosCounter + tMatsToReceive2( Ik ).length() -1 }, { 0, 0 } ) = tMatsToReceive2( Ik ).data();
                    aBbbb ( {tAdofPosCounter, tAdofPosCounter +  tSharedAdofPosLocal( Ik ).length() -1 }, { 0, 0 } ) = tSharedAdofPosLocal( Ik ).data();

                    tAdofPosCounter =tAdofPosCounter + tMatsToReceive2( Ik ).length();
                }
            }
        }
    };

    //-----------------------------------------------------------------------------------------------------------------

    void communicate_check_if_owned_adof_exists( moris::Cell< moris::Cell < Adof * > > & tAdofListofTypes )
    {
        moris::Mat< moris::uint > aCommunicationList( par_size(), 1 );

        for ( moris::uint Ik = 0; Ik < aCommunicationList.length(); Ik++ )
        {
            aCommunicationList( Ik, 0 ) = Ik;
        }

        // Loop over all different adof types in this temporary list
        for ( moris::uint Ij = 0; Ij < tAdofListofTypes.size(); Ij++ )
        {
            moris::Cell< moris::Mat< moris::uint > > tSharedAdofPosGlobal( par_size() );    //What to ask

            // Set Mat to store number of shared adofs per processor
            moris::Mat< moris::uint > tNumSharedAdofsPerProc( par_size(), 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ib = 0; Ib < tAdofListofTypes( Ij ).size(); Ib++ )
            {
                // Check if adof at this position is not NULL
                if ( tAdofListofTypes( Ij )( Ib ) != NULL )
                {
                    // Check if owning processor is this processor
                    if (  tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor() != par_rank() )
                    {
                        // get owning procssor
                        moris::uint tProcID = tAdofListofTypes( Ij )( Ib )->get_adof_owning_processor();

                        // Add +1 to the processor number of shared dofs per processor
                        tNumSharedAdofsPerProc( tProcID, 0) = tNumSharedAdofsPerProc( tProcID, 0) + 1;
                    }
                }
            }

            // Set size of the moris::Mats in the Cell
            for ( moris::uint Ik = 0; Ik < par_size(); Ik++ )
            {
                if ( tNumSharedAdofsPerProc( Ik, 0 ) != 0 )
                {
                    tSharedAdofPosGlobal( Ik ).set_size( tNumSharedAdofsPerProc( Ik, 0 ), 1);
                }
            }

            // Temporary Mat to add external adof ids at the next spot in the matrix which will be communicated
            moris::Mat< moris::uint > tShredAdofPosPerProc( par_size(), 1, 0 );

            // Loop over adofs per type
            for ( moris::uint Ia = 0; Ia < tAdofListofTypes( Ij ).size(); Ia++ )
            {
                // Check if adof at this position is not NULL
                if ( tAdofListofTypes( Ij )( Ia ) != NULL )
                {
                    // Check if owning processor is this processor
                    if (  tAdofListofTypes( Ij )( Ia )->get_adof_owning_processor() != par_rank() )
                    {
                        moris::uint tProcId = tAdofListofTypes( Ij )( Ia )->get_adof_owning_processor();

                        // Add owning procesor id to moris::Mat
                        tSharedAdofPosGlobal( tProcId )( tShredAdofPosPerProc( tProcId, 0 ), 0 ) = tAdofListofTypes( Ij )( Ia )->get_adof_external_id();

                        tShredAdofPosPerProc( tProcId, 0 ) = tShredAdofPosPerProc( tProcId, 0 ) + 1;
                    }
                }
            }

            // receiving list
            moris::Cell< moris::Mat< moris::uint > > tMatsToReceive;

            barrier();

            communicate_mats( aCommunicationList,
                              tSharedAdofPosGlobal,
                              tMatsToReceive );

            // Loop over all Mats set dummy owned adofs
            for ( moris::uint Ik = 0; Ik < tMatsToReceive.size(); Ik++ )
            {
                for ( moris::uint Ii = 0; Ii < tMatsToReceive( Ik ).length(); Ii++ )
                {
                    if ( tAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) == NULL )
                    {
                        tAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) ) = new Adof();
                        tAdofListofTypes( Ij )( tMatsToReceive( Ik )( Ii ) )->set_adof_owning_processor( par_rank() );
                    }
                }
            }
        }
    };

    //-----------------------------------------------------------------------------------------------------------
    void set_pdof_t_matrix()
    {
        // Get number of pdof hosts in pdof host list
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Tell pdofs to get t their t matrix
        for ( moris::uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
        {
            // all pdofs of this pdof host will ask for their t matrix
            mPdofHostList(Ij)->set_t_matrix();
        }
    };

    //-----------------------------------------------------------------------------------------------------------

    const moris::uint get_num_adofs()
    {
        moris::uint tNumAdofs = mAdofList.size();
        return tNumAdofs;
    };

    const moris::Mat< int > get_local_adof_ids()
    {
        Mat< int > tLocalAdofIds ( mAdofList.size(), 1 );
        for ( moris::uint Ij = 0; Ij < mAdofList.size(); Ij++ )
        {
            tLocalAdofIds( Ij, 0 ) = mAdofList( Ij )->get_adof_id();
        }
        return tLocalAdofIds;
    };
};
}
}

#endif /* SRC_FEM_CL_DOF_MANAGER_HPP_ */
