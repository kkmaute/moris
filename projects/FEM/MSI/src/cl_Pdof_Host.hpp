/*
 * cl_Pdof_Host.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_PDOF_HOST_HPP_
#define SRC_FEM_CL_PDOF_HOST_HPP_

#include "cl_Dof_Type_Enums.hpp"
#include "cl_Node_Obj.hpp"
#include "cl_Adof.hpp"

#include "fn_unique.hpp" // LNA/src
#include "cl_Map.hpp" // LNA/src

namespace moris
{
    namespace MSI
    {
    struct Pdof
    {
        moris::uint                  mDofTypeIndex;
        moris::uint                  mTimeStepIndex;
        moris::Mat < moris::sint >   mAdofIds;
        moris::Mat < moris::real >   mTmatrix;

        moris::Cell < Adof* >        mAdofPtrList;
    };

//-------------------------------------------------------------------------------------------------
    class Pdof_Host
    {
    private:
        moris::Cell< enum Dof_Type >        mPdofTypeList;
        moris::Cell< moris::Cell< Pdof* > > mListOfPdofTypeTimeLists;          // FIXME better name and perhaps Cell of Cell

        moris::Mat< moris::uint >               mUniqueAdofList;

        moris::map < moris::uint, moris::uint > mUniqueAdofMap;

    protected:
        mtk::Vertex *  mNodeObj;
        moris::luint mNodeID;

       //FIXME Add interpolation order

    public:
        Pdof_Host()
        {
        };

        Pdof_Host( mtk::Vertex * aNodeObj ) : mNodeObj( aNodeObj )
        {
             mNodeID = mNodeObj->get_id();
        };

        ~Pdof_Host()
        {};

        mtk::Vertex * const get_node_obj_ptr()
        {
            return mNodeObj;
        };

//-------------------------------------------------------------------------------------------------
        //FIXME rewrite this member function
        void set_pdof_type( const enum Dof_Type                  aDof_Type,
                            const moris::Mat< moris::uint >    & aTimeSteps,
                                  moris::Cell< enum Dof_Type > & aPdofTypeList)
        {
            // set a pdof type which belongs to this pdof host.
            bool tDofTypeExists = false;
            bool tDofTypeExistsGlobal = false;
            moris::sint tDofTypeIndex = -1;
            moris::sint tDofTypeIndexGlobal = -1;

            // check if dof type exists. if it exists set flag to true
            for ( moris::uint Ii = 0; Ii < mPdofTypeList.size(); Ii++ )
            {
                if ( mPdofTypeList( Ii ) == aDof_Type )
                {
                    // Check if dof tupe exists twice //FIXMe might be useles after introducing the break
                    MORIS_ERROR( ! tDofTypeExists, "Pdof_Host::set_dof_type(): Dof type exists twice. This is not allowed");

                    tDofTypeExists = true;
                    tDofTypeIndex = Ii;
                    break;
                }
            }

            // check if dof type exists. if it exists set flag to true
            for ( moris::uint Ik = 0; Ik < aPdofTypeList.size(); Ik++ )
            {
                if ( aPdofTypeList( Ik ) == aDof_Type )
                {
                    tDofTypeExistsGlobal = true;
                    tDofTypeIndexGlobal = Ik;
                    break;
                }
            }
            // If dof type does not exist globally, drop error
            MORIS_ERROR(  tDofTypeExistsGlobal, "Pdof_Host::set_dof_type(): Dof type does not exist globally");

            // if dof type does not exist set new dof type.
            if ( !tDofTypeExists )
            {
                tDofTypeIndex = mPdofTypeList.size();

                mPdofTypeList.resize( tDofTypeIndex + 1);
                mPdofTypeList( tDofTypeIndex ) = aDof_Type;

                // FIXME check if resize works
                mListOfPdofTypeTimeLists.resize( tDofTypeIndex + 1);
                mListOfPdofTypeTimeLists( tDofTypeIndex ).resize( aTimeSteps.length() );

                // FIXME use time steps in moris mat dynamically
                // Create new dof type. Add index and time
                mListOfPdofTypeTimeLists( tDofTypeIndex )( aTimeSteps( 0, 0) ) = new Pdof;
                mListOfPdofTypeTimeLists( tDofTypeIndex )( aTimeSteps( 0, 0) ) -> mDofTypeIndex = tDofTypeIndexGlobal;
                mListOfPdofTypeTimeLists( tDofTypeIndex )( aTimeSteps( 0, 0) ) -> mTimeStepIndex = 0;    //FIXME
            }
            else
            {
                MORIS_ERROR( aTimeSteps.length() == mListOfPdofTypeTimeLists( tDofTypeIndex ).size(), " Pdof_Host::set_pdof_type(). Time Levels are not consistent.");
            }
            //std::cout<<"type_Cell "<< mPdofTypeList( 0 )<<std::endl;

            // FIXME return pointer to pdof
        };

//-------------------------------------------------------------------------------------------------
        void get_adofs( moris::Cell< moris::Cell< Adof * > > & aAdofList )
        {
            //Get number of pdof Types in this pdof host
            moris::uint tNumPdofTypes = mListOfPdofTypeTimeLists.size();

            // Loop over all pdof types to create adofs
            for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Add loop for more timesteps Fixme add integration order
                // Get mesh Ids for the used adofs
                moris::Mat < moris::sint > tAdofMeshIds = mNodeObj->get_adof_ids();  //FIXME FIXME FIXME need more information about time and type

                // since petsc requires int, the owner matrix must be casted

                auto tOwners = mNodeObj->get_adof_owners();

                moris::uint tNumberOfOwners = tOwners.length();
                moris::Mat < moris::sint > tAdofOwningProcessorList( tNumberOfOwners, 1 );
                for( uint k=0; k<tNumberOfOwners; ++k )
                {
                	tAdofOwningProcessorList( k ) = tOwners( k );
                }

                // Set size of vector with adpf ptr
                mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofPtrList.resize( tAdofMeshIds.length() );

                // Get pdof type Index
                moris::uint tPdofTypeIndex = mListOfPdofTypeTimeLists( Ii )( 0 )->mDofTypeIndex;

                // loop over all adofs in the matrix and create an adof if it does not exist, yet.
                for ( moris::uint Ik = 0; Ik < tAdofMeshIds.length(); Ik++ )
                {
                    // Check if adof exists
                    if ( aAdofList( tPdofTypeIndex )( tAdofMeshIds( Ik ) ) == NULL)
                    {
                        //std::cout<<*(aAdofList( tPdofTypeIndex )( tAdofMeshIds( 5 ) ))<<std::endl;
                        aAdofList( tPdofTypeIndex )( tAdofMeshIds( Ik ) ) = new Adof();
                        aAdofList( tPdofTypeIndex )( tAdofMeshIds( Ik ) )->set_adof_owning_processor( tAdofOwningProcessorList( Ik ) );
                    }
                    else
                    { }

                    // set pointer to adof on corresponding pdof/time
                    mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofPtrList( Ik ) = aAdofList( tPdofTypeIndex )( tAdofMeshIds( Ik ) );
                }
            }
        };

//-------------------------------------------------------------------------------------------------
        void get_adofs_ids(  )
        {
            //Get number of pdof Types in this pdof host
            moris::uint tNumPdofTypes = mListOfPdofTypeTimeLists.size();

            // Loop over all pdof types to create adofs
            for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Get number of adofs ptr on this pdof/time
                moris::uint tNumAdofPtr = mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofPtrList.size();

                // Set size of matrix containing this pdof/time adof Ids
                mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofIds.set_size( tNumAdofPtr, 1 );

                // loop over all adof ptr of this pdof/time and add the adof Ids to this pdof
                for ( moris::uint Ik = 0; Ik < tNumAdofPtr; Ik++ )
                {
                    mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofIds( Ik, 0 ) = mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofPtrList( Ik )->get_adof_id();
                }
            }
        };

//-------------------------------------------------------------------------------------------------
        void create_unique_adof_list()
        {
            //Get number of pdof Types in this pdof host
            moris::uint tNumPdofTypes = mListOfPdofTypeTimeLists.size();

            moris::uint tAdofCounter = 0;
            // Loop over all adofs of this pdof host to determine maximal number of adofs
            for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++)
            {
                tAdofCounter =tAdofCounter + mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofIds.length();
            }

            moris::Mat< moris::uint > tUniqueAdofList( tAdofCounter, 1 );

            moris::uint tCounter = 0;

            // Loop over all adofs of this pdof host and create a list
            for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++)
            {
                for ( moris::uint Ik = 0; Ik < mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofIds.length(); Ik++)
                {
                    tUniqueAdofList( tCounter, 0 ) = mListOfPdofTypeTimeLists( Ii )( 0 )->mAdofIds( Ik, 0 );
                    tCounter = tCounter + 1;
                }
            }

            // make list unique
            mUniqueAdofList = moris::unique( tUniqueAdofList );
        };
 //-------------------------------------------------------------------------------------------------
         void set_t_matrix()
         {
             //Get number of pdof Types in this pdof host
             moris::uint tNumPdofTypes = mListOfPdofTypeTimeLists.size();

             // Loop over all pdof types to add T matrices
             for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
             {
                 // Add loop for more timesteps
                 auto tTmatrix = mNodeObj->get_t_matrix();
                 mListOfPdofTypeTimeLists( Ii )( 0 )->mTmatrix = tTmatrix->data();
             }
         };

//-------------------------------------------------------------------------------------------------
         void set_unique_adof_map()
         {
             //Get number of unique adofs of this equation object
             moris::uint tNumUniqueAdofs = mUniqueAdofList.length();

             // Loop over all unique adofs of this equation object
             for ( moris::uint Ii = 0; Ii < tNumUniqueAdofs; Ii++ )
             {
                 mUniqueAdofMap[ mUniqueAdofList( Ii, 0 ) ] = Ii;
             }
         };

//-------------------------------------------------------------------------------------------------
        const moris::uint get_num_pdofs()
        {
            //Get number of pdof Types in this pdof host
            moris::uint tNumPdofTypes = mListOfPdofTypeTimeLists.size();

            moris::uint counter = 0;
            // Loop over all pdof types
            for ( moris::uint Ii = 0; Ii < tNumPdofTypes; Ii++ )
            {
                // Add loop for more timesteps
                counter = counter + mListOfPdofTypeTimeLists( Ii ).size();
            }

            return counter;
        };

//-------------------------------------------------------------------------------------------------
        moris::Cell< moris::Cell< Pdof* > > & get_pdof_hosts_pdof_list()
        {
            return mListOfPdofTypeTimeLists;
        }

//-------------------------------------------------------------------------------------------------
        void set_pointer_to_Tmatrix()
        {};

//-------------------------------------------------------------------------------------------------
        void set_adof_IDs()
        {};

//-------------------------------------------------------------------------------------------------
    };


    }
}



#endif /* SRC_FEM_CL_PDOF_HOST_HPP_ */
