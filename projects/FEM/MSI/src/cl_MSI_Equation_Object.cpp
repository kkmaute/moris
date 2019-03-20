/*
 * cl_Equation_Object.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"

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
        mTimeSteps( 0, 0 ) = 1;
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
                                                const Matrix< DDSMat >           & aPdofTypeMap,
                                                      moris::Cell< Pdof_Host * > & aPdofHostList )
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
                mMyPdofHosts( Ii )->set_pdof_type( mEqnObjDofTypeList( Ik ), mTimeSteps( Ik ), aNumUsedDofTypes, aPdofTypeMap );
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

        // Loop over all pdof hosts and dof types. Appending the pdof pointers to the pdof list of this equation object
        for ( moris::uint Ik = 0; Ik < tNumMyPdofHosts; Ik++ )
        {
            // Loop over all pdof types
            for ( moris::uint Ij = 0; Ij < ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() ).size(); Ij++ )
            {
                // Append all time levels of this pdof type
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

         MORIS_ASSERT( tNumUniqueAdofs != 0,"Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. ");

         // Get MAX number of pdofs for this equation object
         moris::uint tNumMyPdofs = mFreePdofs.size();

         MORIS_ASSERT( tNumMyPdofs != 0,"Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. ");


         aPADofMap.set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

         // Loop over all pdofs of this equation object
         for ( moris::uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
         {
             auto tPdof = mFreePdofs( Ii );

             // Loop over all adof Ids of this pdof
             for ( moris::uint Ik = 0; Ik < tPdof->mAdofIds.length(); Ik++ )
             {
                 // Getting tPADofMap column entry for the corresponding value
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

        //print(trans( tTMatrix ), "trans( tTMatrix )");
        //print(mResidual, "mResidual");
        aEqnObjRHS = trans( tTMatrix ) * mResidual;
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::get_my_pdof_values( )
    {
        Matrix< DDRMat > tTMatrix;

        // build T-matrix
        this->build_PADofMap( tTMatrix );

        Matrix< DDRMat > tMyValues;

        // Extract this equation objects adof values from solution vector
        mSolVec->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

        // multiply t_matrix with adof values to get pdof values
        mPdofValues = tTMatrix * tMyValues;

        this->set_vector_entry_number_of_pdof();
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::get_my_pdof_values( const moris::Cell< enum Dof_Type > & aRequestedDofTypes,
                                                    Matrix< DDRMat >             & aRequestedPdofValues )
    {
        // Initialize list which contains the maximal number of time levels per dof type
        Matrix< DDSMat > tTimeLevelsPerDofType( aRequestedDofTypes.size(), 1, -1 );

        moris::sint tCounter = 0;

        // Loop over requested dof types
        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        {
            // Loop over all elemental pdof hosts
            for ( moris::uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
            {
                // Get dof type index
                moris::sint tDofTypeIndex = mModelSolverInterface->get_dof_manager()
                                                                 ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                MORIS_ERROR( mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) !=0,
                        "Equation_Object::get_my_pdof_values: talk with Mathias about this");                         //FIXME delete this error after a closer look

                // get number of time levels for this dof type
                moris::sint tNumTimeLevels = mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex );
                tCounter = tCounter + tNumTimeLevels;

                // Add maximal value of time levels to list
                tTimeLevelsPerDofType( Ii, 0 ) = std::max( tTimeLevelsPerDofType( Ii, 0 ), tNumTimeLevels );
            }
            MORIS_ASSERT( tTimeLevelsPerDofType( Ii, 0 ) > -1, "Equation_Object::get_my_pdof_values: no time levels exist on this dof type on element %-5i", mEqnObjInd );
        }
        // Set size matrix for requestes pdof values
        aRequestedPdofValues.resize( tCounter, 1 );

        moris::sint tCounter_2 = 0;

        // Loop over requested dof types
        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        {
            // Get maximal Number of time levels on this pdof type
            moris::sint tMaxTimeLevelsOnDofType = tTimeLevelsPerDofType( Ii, 0 );

            // Loop over this pdofs time levels
            for ( moris::sint Ia = 0; Ia < tMaxTimeLevelsOnDofType; Ia++ )
            {
                // Loop over all elemental pdof hosts
                for ( moris::uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
                {
                    // Get dof type index
                    moris::sint tDofTypeIndex = mModelSolverInterface->get_dof_manager()
                                                                     ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                    // Check if number if time levels on this dof type is smaller than maximal number of time levels on dof type
                    if ( (sint)mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) == tMaxTimeLevelsOnDofType )
                    {
                        // get pointer list all time pdofs on this pdof type
                        moris::Cell< Pdof* > tPdofTimeList = mMyPdofHosts( Ik )->get_pdof_time_list( tDofTypeIndex );

                        // get entry number of this pdof in the elemental pdof value vector
                        moris::uint mElementalSolVecEntry = tPdofTimeList( Ia )->mElementalSolVecEntry;

                        // Put this pdof value into the requested pdof vector
                        aRequestedPdofValues( tCounter_2++, 0 ) = mPdofValues( mElementalSolVecEntry , 0 );
                    }
                }
            }
        }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::set_vector_entry_number_of_pdof()
    {
        moris::uint tNumMyPdofs = mFreePdofs.size();
        // Loop over all pdofs of this element
        for ( moris::uint Ik = 0; Ik < tNumMyPdofs; Ik++ )
        {
            mFreePdofs( Ik )->mElementalSolVecEntry = Ik;
        }
    }
}
}
