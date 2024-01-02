/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Pdof_Host.cpp
 *
 */

#include "cl_MSI_Pdof_Host.hpp"
#include "cl_FEM_Node_Base.hpp"

#include "cl_MSI_Model_Solver_Interface.hpp"

#include "fn_isfinite.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace MSI
    {
        //-----------------------------------------------------------------------------------------------------------

        Pdof_Host::Pdof_Host(
                const uint      aNumUsedDofTypes,
                fem::Node_Base* aNodeObj )
                : mNodeObj( aNodeObj )
        {
            mNodeID = mNodeObj->get_id();
            // mNodeInd = mNodeObj->get_index();

            mPdofTypeExist.set_size( aNumUsedDofTypes, 1, 0 );

            // Set size of list to the number of used nodes
            mListOfPdofTimePerType.resize( aNumUsedDofTypes );
        }

        //-----------------------------------------------------------------------------------------------------------

        Pdof_Host::~Pdof_Host()
        {
            for ( Vector< Pdof* >& tList : mListOfPdofTimePerType )
            {
                for ( Pdof* tPdof : tList )
                {
                    delete tPdof;
                }
                tList.clear();
            }
            mListOfPdofTimePerType.clear();
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::set_pdof_type(
                const enum Dof_Type     aDof_Type,
                const Matrix< DDUMat >& aTimePerDofType,
                const uint              aNumUsedDofTypes,
                const Matrix< DDSMat >& aPdofTypeMap )
        {
            // Get global dof type index
            sint tDofTypeIndex = aPdofTypeMap( static_cast< int >( aDof_Type ) );

            // Get number of time levels for dof type
            uint tNumTimeLevelforDofType = aTimePerDofType( tDofTypeIndex, 0 );

            // if dof type does not exist set new dof type.
            if ( mPdofTypeExist( tDofTypeIndex ) == 0 )
            {
                // Set mPdofTypeExist to 1. ==> Dof type exists
                mPdofTypeExist( tDofTypeIndex ) = 1;

                mListOfPdofTimePerType( tDofTypeIndex ).resize( tNumTimeLevelforDofType, nullptr );

                // Create new dof type. Add index and time
                for ( uint Ii = 0; Ii < tNumTimeLevelforDofType; Ii++ )
                {
                    // Create pdof
                    mListOfPdofTimePerType( tDofTypeIndex )( Ii ) = new Pdof;

                    // Set dof type index
                    mListOfPdofTimePerType( tDofTypeIndex )( Ii )->mDofTypeIndex = tDofTypeIndex;

                    //
                    mListOfPdofTimePerType( tDofTypeIndex )( Ii )->mTimeStepIndex = Ii;
                }
            }
            else
            {
                MORIS_ASSERT( tNumTimeLevelforDofType == mListOfPdofTimePerType( tDofTypeIndex ).size(), " Pdof_Host::set_pdof_type(). Time Levels are not consistent." );
            }

            // FIXME return pointer to pdof?
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::get_adofs(
                const Matrix< DDUMat >&              aTimeLevelOffsets,
                Vector< Vector< Adof* > >& aAdofList,
                Model_Solver_Interface*              aModelSolverInterface,
                const bool&                          aUseHMR )
        {
            if ( aUseHMR )
            {
                this->create_adofs_based_on_Tmatrix( aTimeLevelOffsets, aAdofList, aModelSolverInterface );
            }
            else
            {
                this->create_adofs_based_on_pdofs( aTimeLevelOffsets, aAdofList );
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::create_adofs_based_on_Tmatrix(
                const Matrix< DDUMat >&              aTimeLevelOffsets,
                Vector< Vector< Adof* > >& aAdofList,
                Model_Solver_Interface*              aModelSolverInterface )
        {
            // Get number of DoF Types (i.e. number of PDofs created on this node/host per DoF in time direction)
            uint tNumPdofTypes = mListOfPdofTimePerType.size();

            // Loop over all pdof types to create adofs
            for ( uint iPDofType = 0; iPDofType < tNumPdofTypes; iPDofType++ )
            {
                // Ask for adof mesh index for this dof type
                uint tAdofMeshIndex = (uint)aModelSolverInterface->get_adof_index_for_type( iPDofType );

                if ( mListOfPdofTimePerType( iPDofType ).size() != 0 )
                {
                    // Get mesh Ids for the used adofs
                    Matrix< DDSMat > tAdofMeshId  = mNodeObj->get_adof_ids( tAdofMeshIndex );
                    Matrix< DDSMat > tAdofMeshInd = mNodeObj->get_adof_indices( tAdofMeshIndex );
                    Matrix< IdMat >  tOwners      = mNodeObj->get_adof_owners( tAdofMeshIndex );

                    // check that the node (i.e. the PDoF host) actually has an interpolation rule
                    MORIS_ASSERT(
                            tAdofMeshInd.numel() != 0,
                            "Pdof_Host::create_adofs_based_on_Tmatrix() - "
                            "Node with T-matrix size of 0 at coords %s which is owned by proc #%i.",
                            ios::stringify_log( mNodeObj->get_geometric_vertex()->get_coords() ).c_str(),
                            mNodeObj->get_geometric_vertex()->get_owner() );

                    for ( uint iTimeLevel = 0; iTimeLevel < mListOfPdofTimePerType( iPDofType ).size(); iTimeLevel++ )
                    {
                        // Set size of vector with adof pointer
                        mListOfPdofTimePerType( iPDofType )( iTimeLevel )->mAdofPtrList.resize( tAdofMeshInd.numel() );

                        // Get pdof type Index
                        uint tPdofTypeIndex = mListOfPdofTimePerType( iPDofType )( iTimeLevel )->mDofTypeIndex;

                        uint tAdofType = aTimeLevelOffsets( tPdofTypeIndex );

                        // loop over all adofs in the matrix and create an adof if it does not exist, yet.
                        for ( uint iIpMesh = 0; iIpMesh < tAdofMeshInd.numel(); iIpMesh++ )
                        {
                            sint tCurrentADofMeshIndex = tAdofMeshInd( iIpMesh );
                            uint tADofTypeWithTime     = tAdofType + iTimeLevel;

                            // if adof does not exist, store it based on adof type & time and generic coefficient
                            // index (here: tCurrentADofMeshIndex) for particular mesh
                            if ( aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex ) == nullptr )
                            {
                                // Create new adof pointer. Put adof on the right spot of the temporary vector
                                aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex ) = new Adof();

                                // Set this adofs owning processor
                                aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex )->set_adof_owning_processor( tOwners( iIpMesh ) );

                                // Set adof external Id and Ind. Id used for comm, Ind used for HMR ordering
                                aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex )->set_adof_external_id( tAdofMeshId( iIpMesh ) );    // FIXME delete
                                aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex )->set_adof_external_ind( tAdofMeshInd( iIpMesh ) );
                            }

                            // set pointer to adof on corresponding pdof/time
                            mListOfPdofTimePerType( iPDofType )( iTimeLevel )->mAdofPtrList( iIpMesh ) = aAdofList( tADofTypeWithTime )( tCurrentADofMeshIndex );
                        }
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::create_adofs_based_on_pdofs(
                const Matrix< DDUMat >&              aTimeLevelOffsets,
                Vector< Vector< Adof* > >& aAdofList )
        {
            //        //Get number of pdof Types in this pdof host
            //        uint tNumPdofTypes = mListOfPdofTimePerType.size();
            //
            //        // Loop over all pdof types to create adofs
            //        for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            //        {
            //            if ( mListOfPdofTimePerType( Ii ).size() != 0 )
            //            {
            //                 // Get mesh Ids for the used adofs
            //                 sint tAdofMeshId = mNodeObj->get_id();
            //                 sint tAdofMeshInd = mNodeObj->get_index();
            //
            //                 // since petsc requires int, the owner matrix must be casted
            //                 auto tOwner = mNodeObj->get_owner();
            //
            //                 for ( uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
            //                 {
            //                    // Set size of vector with adpf ptr
            //                    mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList.resize( 1 );
            //
            //                    // Get pdof type Index
            //                    uint tPdofTypeIndex = mListOfPdofTimePerType( Ii )( Ij )->mDofTypeIndex;                  ///////
            //
            //                    uint tAdofType = aTimeLevelOffsets( tPdofTypeIndex, 0 );
            //
            //                    // Check if adof exists
            //                    if ( aAdofList( tAdofType + Ij )( tAdofMeshInd ) == NULL)
            //                    {
            //                        // Create new adof pointer. Put adof on the right spot of the temporary vector
            //                        aAdofList( tAdofType + Ij )( tAdofMeshInd ) = new Adof();
            //
            //                        // Set this adofs owning processor
            //                        aAdofList( tAdofType + Ij )( tAdofMeshInd )->set_adof_owning_processor( tOwner );
            //
            //                        // Set adof external Id and Ind. Id used for comm, Ind used for HMR ordering
            //                        aAdofList( tAdofType + Ij )( tAdofMeshInd )->set_adof_external_id( tAdofMeshId );               //FIXME delete
            //                        aAdofList( tAdofType + Ij )( tAdofMeshInd )->set_adof_external_ind( tAdofMeshInd );
            //                    }
            //
            //                    // set pointer to adof on corresponding pdof/time
            //                    mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList( 1 ) = aAdofList( tAdofType + Ij )( tAdofMeshInd );
            //                }
            //            }
            //        }
        }

        //-----------------------------------------------------------------------------------------------------------
        void
        Pdof_Host::get_adofs_ids()
        {
            // Get number of pdof Types in this pdof host
            uint tNumPdofTypes = mListOfPdofTimePerType.size();

            // Loop over all pdof types
            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Loop over all timelevel
                for ( uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
                {
                    // Get number of adofs ptr on this pdof/time
                    uint tNumAdofPtr = mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList.size();

                    // Set size of matrix containing this pdof/time adof Ids
                    mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.set_size( tNumAdofPtr, 1 );

                    // loop over all adof ptr of this pdof/time and add the adof Ids to this pdof
                    for ( uint Ik = 0; Ik < tNumAdofPtr; Ik++ )
                    {
                        auto tPointer                                         = mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList( Ik );
                        mListOfPdofTimePerType( Ii )( Ij )->mAdofIds( Ik, 0 ) = tPointer->get_adof_id();
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        //        void Pdof_Host::create_unique_adof_list()
        //        {
        //            //Get number of pdof Types in this pdof host
        //            uint tNumPdofTypes = mListOfPdofTimePerType.size();
        //
        //            uint tAdofCounter = 0;
        //            // Loop over all adofs of this pdof host to determine maximal number of adofs
        //            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++)
        //            {
        //                for ( uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
        //                {
        //                    tAdofCounter = tAdofCounter + mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.numel();
        //                }
        //            }
        //
        //            Matrix< DDUMat > tUniqueAdofList( tAdofCounter, 1 );
        //
        //            uint tCounter = 0;
        //
        //            // Loop over all adofs of this pdof host and create a list of adof ids
        //            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++)
        //            {
        //                for ( uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
        //                {
        //                    for ( uint Ik = 0; Ik < mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.numel(); Ik++)
        //                    {
        //                        tUniqueAdofList( tCounter, 0 ) = mListOfPdofTimePerType( Ii )( Ij )->mAdofIds( Ik, 0 );
        //                        tCounter++;
        //                    }
        //                }
        //            }
        //
        //            // make list unique
        //            moris::unique( tUniqueAdofList, mUniqueAdofList );
        //        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::set_t_matrix(
                const bool&             aUseHMR,
                Model_Solver_Interface* aModelSolverInterface )
        {
            // Get number of pdof Types in this pdof host
            uint tNumPdofTypes = mListOfPdofTimePerType.size();

            // Loop over all pdof types and times to add T matrices
            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Ask for adof mesh index for this dof type
                uint tAdofMeshIndex = (uint)aModelSolverInterface->get_adof_index_for_type( Ii );

                for ( uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
                {
                    if ( aUseHMR )
                    {
                        // Get TMatrix. Add Tmatrix to type and time list
                        const Matrix< DDRMat >* tTmatrix             = mNodeObj->get_t_matrix( tAdofMeshIndex );
                        mListOfPdofTimePerType( Ii )( Ij )->mTmatrix = tTmatrix->matrix_data();
                    }
                    else
                    {
                        mListOfPdofTimePerType( Ii )( Ij )->mTmatrix.resize( 1, 1 );
                        mListOfPdofTimePerType( Ii )( Ij )->mTmatrix( 0, 0 ) = 1.0;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        uint
        Pdof_Host::get_num_pdofs()
        {
            // Get number of pdof Types in this pdof host
            uint tNumPdofTypes = mListOfPdofTimePerType.size();

            uint counter = 0;
            // Loop over all pdof types
            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Add loop for more time steps
                counter = counter + mListOfPdofTimePerType( Ii ).size();
            }
            return counter;
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Pdof_Host::print_t_matrix(
                const bool&             aUseHMR,
                Model_Solver_Interface* aModelSolverInterface )
        {
            // Get number of pdof Types in this pdof host
            uint tNumPdofTypes = mListOfPdofTimePerType.size();

            // Array for coordinates
            Matrix< DDRMat > tCoords;

            // Loop over all pdof types and times to add T matrices
            for ( uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Ask for adof mesh index for this dof type
                uint tAdofMeshIndex = (uint)aModelSolverInterface->get_adof_index_for_type( 0 );

                // nothing is printed unless Tmatrix is provided by HMR
                if ( aUseHMR )
                {
                    // Get TMatrix.
                    const Matrix< DDRMat >* tTmatrix    = mNodeObj->get_t_matrix( tAdofMeshIndex );
                    const Matrix< DDSMat >  tAdofMeshId = mNodeObj->get_adof_ids( tAdofMeshIndex );

                    // Get nodal coordinates
                    mNodeObj->get_vertex_coords( tCoords );

                    std::cout << "%=======================================================================================\n";

                    std::cout << "% Adof_index: " << Ii << std::endl;

                    print( tCoords, "Pdof_Host_Coordinates" );

                    print( *tTmatrix, "Tmatrix weights" );

                    print( tAdofMeshId, "Tmatrix Ids" );

                    std::cout << "%=======================================================================================\n";
                }
            }
        }
    }    // namespace MSI
}    // namespace moris
