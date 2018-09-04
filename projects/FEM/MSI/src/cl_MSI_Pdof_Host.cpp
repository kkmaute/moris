/*
 * cl_Pdof_Host.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_FEM_Node_Base.hpp"

namespace moris
{
namespace MSI
{
    Pdof_Host::Pdof_Host( const moris::uint      aNumUsedDofTypes,
                                fem::Node_Base * aNodeObj ) : mNodeObj( aNodeObj )
    {
        mNodeID = mNodeObj->get_id();
        //mNodeInd = mNodeObj->get_index();

        mPdofTypeExist.set_size( aNumUsedDofTypes, 1, 0 );

        // Set size of list to the number of used nodes
        mListOfPdofTimePerType.resize( aNumUsedDofTypes );
    }

    Pdof_Host::~Pdof_Host()
    {
        for ( moris::uint Ik = 0; Ik < mListOfPdofTimePerType.size(); Ik++ )
        {
            for ( moris::uint Ii = 0; Ii < mListOfPdofTimePerType( Ik ).size(); Ii++ )
            {
                delete mListOfPdofTimePerType( Ik )( Ii );
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Pdof_Host::set_pdof_type( const enum Dof_Type                  aDof_Type,
                                   const moris::Mat< moris::uint >    & aTimeSteps,
                                   const moris::uint                    aNumUsedDofTypes,
                                   const moris::Mat< moris::sint >    & aPdofTypeMap)
    {
        // Get global dof type index
        moris::sint tDofTypeIndex = aPdofTypeMap( static_cast< int >( aDof_Type ) );

        // if dof type does not exist set new dof type.
        if( mPdofTypeExist( tDofTypeIndex ) == 0 )
        {
            // Set mPdofTypeExist to 1. ==> Dof type exists
            mPdofTypeExist( tDofTypeIndex ) = 1;

            mListOfPdofTimePerType( tDofTypeIndex ).resize( aTimeSteps.length() );

            // Create new dof type. Add index and time
            for ( moris::uint Ii = 0; Ii < aTimeSteps.length(); Ii++ )
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
            MORIS_ERROR( aTimeSteps.length() == mListOfPdofTimePerType( tDofTypeIndex ).size(), " Pdof_Host::set_pdof_type(). Time Levels are not consistent.");
        }

        // FIXME return pointer to pdof?
    }

    //-----------------------------------------------------------------------------------------------------------
    void Pdof_Host::get_adofs( const moris::Mat< moris::uint >            & aTimeLevelOffsets,
                                     moris::Cell< moris::Cell< Adof * > > & aAdofList )
    {
        //Get number of pdof Types in this pdof host
        moris::uint tNumPdofTypes = mListOfPdofTimePerType.size();

        // Loop over all pdof types to create adofs
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
        {
            if ( mListOfPdofTimePerType( Ii ).size() != 0 )
            {
                 // Get mesh Ids for the used adofs
                 moris::Mat < moris::sint > tAdofMeshId = mNodeObj->get_adof_ids();                      //FIXME add interpolation order in ()
                 moris::Mat < moris::sint > tAdofMeshInd = mNodeObj->get_adof_ind();                      //FIXME add interpolation order in ()

                 // since petsc requires int, the owner matrix must be casted
                 auto tOwners = mNodeObj->get_adof_owners();

                 moris::uint tNumberOfOwners = tOwners.length();

                 moris::Mat < moris::sint > tAdofOwningProcessorList( tNumberOfOwners, 1 );

                 for( uint k=0; k<tNumberOfOwners; ++k )
                 {
                     tAdofOwningProcessorList( k ) = tOwners( k );
                 }

                 for ( moris::uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
                 {
                    // Set size of vector with adpf ptr
                    mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList.resize( tAdofMeshInd.length() );

                    // Get pdof type Index
                    moris::uint tPdofTypeIndex = mListOfPdofTimePerType( Ii )( Ij )->mDofTypeIndex;                  ///////

                    moris::uint tAdofType = aTimeLevelOffsets( tPdofTypeIndex, 0 );

                    // loop over all adofs in the matrix and create an adof if it does not exist, yet.
                    for ( moris::uint Ik = 0; Ik < tAdofMeshInd.length(); Ik++ )
                    {
                        // Check if adof exists
                        if ( aAdofList( tAdofType + Ij )( tAdofMeshInd( Ik ) ) == NULL)
                        {
                            // Create new adof pointer. Put adof on the right spot of the temporary vector
                            aAdofList( tAdofType + Ij )( tAdofMeshInd( Ik ) ) = new Adof();

                            // Set this adofs owning processor
                            aAdofList( tAdofType + Ij )( tAdofMeshInd( Ik ) )->set_adof_owning_processor( tAdofOwningProcessorList( Ik ) );

                            aAdofList( tAdofType + Ij )( tAdofMeshInd( Ik ) )->set_adof_external_id( tAdofMeshId( Ik ) );               //FIXME delete
                        }

                        // set pointer to adof on corresponding pdof/time
                        mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList( Ik ) = aAdofList( tAdofType + Ij )( tAdofMeshInd( Ik ) );
                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Pdof_Host::get_adofs_ids()
    {
        //Get number of pdof Types in this pdof host
        moris::uint tNumPdofTypes = mListOfPdofTimePerType.size();

        // Loop over all pdof types
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
        {
            // Loop over all timelevel
            for ( moris::uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
            {
                // Get number of adofs ptr on this pdof/time
                moris::uint tNumAdofPtr = mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList.size();

                // Set size of matrix containing this pdof/time adof Ids
                mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.set_size( tNumAdofPtr, 1 );

                // loop over all adof ptr of this pdof/time and add the adof Ids to this pdof
                for ( moris::uint Ik = 0; Ik < tNumAdofPtr; Ik++ )
                {
                    auto tPointer = mListOfPdofTimePerType( Ii )( Ij )->mAdofPtrList( Ik );
                    mListOfPdofTimePerType( Ii )( Ij )->mAdofIds( Ik, 0 ) = tPointer->get_adof_id();
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    void Pdof_Host::create_unique_adof_list()
    {
        //Get number of pdof Types in this pdof host
        moris::uint tNumPdofTypes = mListOfPdofTimePerType.size();

        moris::uint tAdofCounter = 0;
        // Loop over all adofs of this pdof host to determine maximal number of adofs
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++)
        {
            for ( moris::uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
            {
                tAdofCounter = tAdofCounter + mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.length();
            }
        }

        moris::Mat< moris::uint > tUniqueAdofList( tAdofCounter, 1 );

        moris::uint tCounter = 0;

        // Loop over all adofs of this pdof host and create a list of adof ids
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++)
        {
            for ( moris::uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
            {
                for ( moris::uint Ik = 0; Ik < mListOfPdofTimePerType( Ii )( Ij )->mAdofIds.length(); Ik++)
                {
                    tUniqueAdofList( tCounter, 0 ) = mListOfPdofTimePerType( Ii )( Ij )->mAdofIds( Ik, 0 );
                    tCounter = tCounter + 1;
                }
            }
        }

        // make list unique
        mUniqueAdofList = moris::unique( tUniqueAdofList );
    }

    //-----------------------------------------------------------------------------------------------------------
    void Pdof_Host::set_t_matrix()
    {
        //Get number of pdof Types in this pdof host
        moris::uint tNumPdofTypes = mListOfPdofTimePerType.size();

        // Loop over all pdof types and times to add T matrices
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
        {
            for ( moris::uint Ij = 0; Ij < mListOfPdofTimePerType( Ii ).size(); Ij++ )
            {
                // Get TMatrix. Add Tmatrix to type and time list
                const moris::Mat< moris::real > * tTmatrix = mNodeObj->get_t_matrix();           //FIXME interpolation order //FIXME FIXME FIXME FIXME FIXME
                mListOfPdofTimePerType( Ii )( Ij )->mTmatrix = tTmatrix->data();
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------
    const moris::uint Pdof_Host::get_num_pdofs()
    {
        //Get number of pdof Types in this pdof host
        moris::uint tNumPdofTypes = mListOfPdofTimePerType.size();

        moris::uint counter = 0;
        // Loop over all pdof types
        for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
        {
            // Add loop for more timesteps
            counter = counter + mListOfPdofTimePerType( Ii ).size();
        }
        return counter;
    }

    //-----------------------------------------------------------------------------------------------------------
}
}
