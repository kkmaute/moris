/*
 * cl_Equation_Object.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Pdof_Host.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    namespace MSI
    {

    Equation_Object::Equation_Object( const moris::Cell< fem::Node_Base * > & aNodeObjs ) : mNodeObj( aNodeObjs )
    {
        mTimeSteps.resize( 1, 1 );
        mTimeSteps( 0, 0 ) = 0;
    }

//-------------------------------------------------------------------------------------------------
    moris::uint Equation_Object::get_max_pdof_hosts_ind()
    {
        auto tMaxPdofHostsInd = mNodeObj( 0 )->get_index();

        // Loop over all node obj. get the maximal node index.
        for ( moris::uint Ii=1; Ii < mNodeObj.size(); Ii++ )
        {
            tMaxPdofHostsInd = std::max( tMaxPdofHostsInd, mNodeObj( Ii )->get_index() );
        }
        return ( moris::uint ) tMaxPdofHostsInd;
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_pdof_hosts( const moris::uint                  aNumUsedDofTypes,
                                                const Matrix< DDSMat >  & aPdofTypeMap,
                                                      moris::Cell< Pdof_Host * > & aPdofHostList)
    {
        // Determine size of list containing this equations objects pdof hosts
        moris::uint tNumMyPdofHosts = mNodeObj.size();                            //Fixme Add ghost and element numbers

        // Resize list containing this equations objects pdof hosts
        mMyPdofHosts.resize( tNumMyPdofHosts, nullptr );

        // Loop over all nodes of this element, creating new pdof hosts if not existing yet.
        for ( moris::uint Ii=0; Ii < mNodeObj.size(); Ii++ )
        {
            // Save node id of node Ii in temporary variable for more clarity.
            //auto tNodeID = mNodeObj( Ii )->get_id();
            auto tNodeID = mNodeObj( Ii )->get_index();

            // check if pdof host corresponding to this node exists.
            if ( aPdofHostList( tNodeID ) == NULL)
            {
                // If node does not exist, create new pdof host.
                aPdofHostList( tNodeID ) = new Pdof_Host( aNumUsedDofTypes, mNodeObj( Ii ) );
            }

            // Add pointer to pdof host to the list containing this equation objects pdof hosts.
            mMyPdofHosts( Ii ) = aPdofHostList( tNodeID );

            // FIXME rewrite this function
            for ( moris::uint Ik=0; Ik < mEqnObjDofTypeList.size(); Ik++ )
            {
                mMyPdofHosts( Ii )->set_pdof_type( mEqnObjDofTypeList( Ik ), mTimeSteps, aNumUsedDofTypes, aPdofTypeMap );
            }
        }

        // Fixme add element
       // FIXME return pointer to pdofs
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_pdof_list()
    {
        // Get number of pdof hosts corresponding to this equation object
        moris::uint tNumMyPdofHosts = mMyPdofHosts.size();

        // Loop over all pdof hosts and get their number of (free) pdofs
        moris::uint tNumMyFreePdofs = 0;
        for ( moris::uint Ii=0; Ii < tNumMyPdofHosts; Ii++ )
        {
            tNumMyFreePdofs = tNumMyFreePdofs + mMyPdofHosts( Ii )->get_num_pdofs();
        }

        // Set size of vector containing this equation objects free pdofs.
        mFreePdofs.reserve( tNumMyFreePdofs );

        // Loop over all pdof hosts and dof types. Appending the psof pointers to the pdof list of this equation object
        for ( moris::uint Ik = 0; Ik < tNumMyPdofHosts; Ik++ )
        {
            // Loop over all pdof types
            for ( moris::uint Ij = 0; Ij < ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() ).size(); Ij++ )
            {
                mFreePdofs.append( ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() )( Ij ) );
            }
        }
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_list_of_adof_ids()
    {
        // Get MAX number of pdofs for this equation object
        moris::uint tNumMyPdofs = mFreePdofs.size();

        // Loop over all pdofs to count their adofs
        moris::uint tNumMyAdofs = 0;
        for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
        {
            // Get Number of adofs cooresponding to this pdof
            moris::uint tNumAdofForThisPdof = ( mFreePdofs( Ij )->mAdofIds ).length();
            tNumMyAdofs = tNumMyAdofs + tNumAdofForThisPdof;
        }

        // Temporary matrix for adofs Ids
        Matrix< DDSMat > tNonUniqueAdofIds( tNumMyAdofs, 1 );

        moris::uint tAdofPosCounter = 0;

        // Loop over all pdofs to get their adofs and put them into a unique list
        for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
        {
            tNonUniqueAdofIds ( {tAdofPosCounter, tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).length() -1 }, { 0, 0} ) = mFreePdofs( Ij )->mAdofIds.matrix_data();

            // Add number if these adofs to number of assembled adofs
            tAdofPosCounter =tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).length();
        }
        // make list of unique Ids
        moris::unique( tNonUniqueAdofIds, mUniqueAdofList );
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::set_unique_adof_map()
    {
        //Get number of unique adofs of this equation object
        moris::uint tNumUniqueAdofs = mUniqueAdofList.length();

        // Loop over all unique adofs of this equation object
        for ( moris::uint Ii = 0; Ii < tNumUniqueAdofs; Ii++ )
        {
            mUniqueAdofMap[ mUniqueAdofList( Ii, 0 ) ] = Ii;
        }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::build_PADofMap( Matrix< DDRMat > & aPADofMap )
    {
         //Get number of unique adofs of this equation object
         moris::uint tNumUniqueAdofs = mUniqueAdofList.length();

         // Get MAX number of pdofs for this equation object
         moris::uint tNumMyPdofs = mFreePdofs.size();

         aPADofMap.set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

         // Loop over all pdofs of this equation object
         for ( moris::uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
         {
             auto tPdof = mFreePdofs( Ii );

             // Loop over all adof Ids of this pdof
             for ( moris::uint Ik = 0; Ik < tPdof->mAdofIds.length(); Ik++ )
             {
                 // Getting tPADofMap column entry for the corrsponding value
                 moris::uint tColumnPos = mUniqueAdofMap[ tPdof->mAdofIds( Ik, 0 ) ];

                 // Insert value into pdof-adof-map
                 aPADofMap( Ii, tColumnPos ) = ( mFreePdofs( Ii )->mTmatrix)( Ik, 0 );
             }
         }
     }

//-------------------------------------------------------------------------------------------------
    moris_index Equation_Object::get_node_index( const moris_index aElementLocalNodeIndex ) const
    {
        return mNodeObj( aElementLocalNodeIndex )->get_index();
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::get_equation_obj_residual( Matrix< DDRMat > & aEqnObjRHS, Dist_Vector * aSolutionVector )
    {
        mSolVec = aSolutionVector;

        this->compute_residual();

        Matrix< DDRMat > tTMatrix;

        this->build_PADofMap( tTMatrix );

        aEqnObjRHS = trans( tTMatrix ) * mResidual;
    }

//-------------------------------------------------------------------------------------------------
}
}
