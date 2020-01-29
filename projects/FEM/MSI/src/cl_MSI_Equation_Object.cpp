/*
 * cl_Equation_Object.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    namespace MSI
    {

    Equation_Object::Equation_Object( const moris::Cell < moris::Cell< fem::Node_Base * > > & aNodeObjs ) : mNodeObj( aNodeObjs )
    {
    }

//-------------------------------------------------------------------------------------------------
    moris::uint Equation_Object::get_max_pdof_hosts_ind()
    {
        auto tMaxPdofHostsInd = mNodeObj( 0 )( 0 )->get_index();

        // Loop over all node obj. get the maximal node index.
        for ( moris::uint Ik = 0; Ik < mNodeObj.size(); Ik++ )
        {
            for ( moris::uint Ii=0; Ii < mNodeObj( Ik ).size(); Ii++ )
            {
                tMaxPdofHostsInd = std::max( tMaxPdofHostsInd, mNodeObj( Ik )( Ii )->get_index() );
            }
        }
        return ( moris::uint ) tMaxPdofHostsInd;
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_pdof_hosts( const moris::uint                  aNumUsedDofTypes,
                                                const Matrix< DDSMat >           & aPdofTypeMap,
                                                const Matrix< DDUMat >           & aTimePerDofType,
                                                      moris::Cell< Pdof_Host * > & aPdofHostList )
    {
//        moris::uint tNumMyPdofHosts = 0;
//
//        // Loop over all node obj. get the maximal node index.
//        for ( moris::uint Ik = 0; Ik < mNodeObj.size(); Ik++ )
//        {
//            // Determine size of list containing this equations objects pdof hosts
//            moris::uint tNumMyPdofHosts = tNumMyPdofHosts + mNodeObj( Ik ).size();              //Fixme Add ghost and element numbers
//        }

        // Resize list containing this equations objects pdof hosts set
        mNumPdofSystems = mNodeObj.size();
        mMyPdofHosts.resize( mNumPdofSystems );                        //Fixme Add ghost and element numbers

        for ( moris::uint Ik = 0; Ik < mNumPdofSystems; Ik++ )
        {
            moris::uint tNumMyPdofHosts = mNodeObj( Ik ).size();

            // Resize list containing this equations objects pdof hosts
            mMyPdofHosts( Ik ).resize( tNumMyPdofHosts, nullptr );

            // Loop over all nodes of this element, creating new pdof hosts if not existing yet.
            for ( moris::uint Ii=0; Ii < mNodeObj( Ik ).size(); Ii++ )
            {
                // Save node id of node Ii in temporary variable for more clarity.
                //auto tNodeID = mNodeObj( Ik )( Ii )->get_id();
                auto tNodeID = mNodeObj( Ik )( Ii )->get_index();

                // check if pdof host corresponding to this node exists.
                if ( aPdofHostList( tNodeID ) == NULL)
                {
                    // If node does not exist, create new pdof host.
                    aPdofHostList( tNodeID ) = new Pdof_Host( aNumUsedDofTypes, mNodeObj( Ik )( Ii ) );
                }

                // Add pointer to pdof host to the list containing this equation objects pdof hosts.
                mMyPdofHosts( Ik )( Ii ) = aPdofHostList( tNodeID );

                // FIXME rewrite this function
                for ( moris::uint Ij=0; Ij < mEquationBlock->get_unique_dof_type_list().size(); Ij++ )
                {
                    mMyPdofHosts( Ik )( Ii )->set_pdof_type( mEquationBlock->get_unique_dof_type_list()( Ij ),
                                                       aTimePerDofType,
                                                       aNumUsedDofTypes,
                                                       aPdofTypeMap );
                }
            }
        }

        // Fixme add element
       // FIXME return pointer to pdofs
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_pdof_list()             //FIXME add time and type
    {
        moris::uint tNumMyFreePdofs = 0;

        // Loop over all pdof hosts and get their number of (free) pdofs
        for ( moris::uint Ik=0; Ik < mMyPdofHosts.size(); Ik++ )
        {
            // Get number of pdof hosts corresponding to this equation object
            moris::uint tNumMyPdofHosts = mMyPdofHosts( Ik ).size();

            for ( moris::uint Ii=0; Ii < tNumMyPdofHosts; Ii++ )
            {
                tNumMyFreePdofs = tNumMyFreePdofs + mMyPdofHosts( Ik )( Ii )->get_num_pdofs();
            }
        }

        // Set size of vector containing this equation objects free pdofs.
        mFreePdofs.reserve( tNumMyFreePdofs );

        for ( moris::uint Ia=0; Ia < mMyPdofHosts.size(); Ia++ )
        {
            moris::uint tNumMyPdofHosts = mMyPdofHosts( Ia ).size();

            // Loop over all pdof types. Ask the first pdof host for the number of pdof types
            for ( moris::uint Ij = 0; Ij < ( mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list() ).size(); Ij++ )
            {
                // Loop over all time levels for this dof type
                for ( moris::uint Ii = 0; Ii < mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list()( Ij ).size(); Ii++ )
                {
                    // Loop over all pdof hosts and dof types. Appending the pdof pointers to the pdof list of this equation object
                    for ( moris::uint Ik = 0; Ik < tNumMyPdofHosts; Ik++ )
                    {
                        // Append all time levels of this pdof type
                        mFreePdofs.push_back( ( mMyPdofHosts( Ia )( Ik )->get_pdof_hosts_pdof_list() )( Ij )( Ii ) );
                    }
                }
            }
        }

        //----------------------------------------------------------------------------------------------------------

        // Ask the first pdof host for the number of pdof types //FIXME
        mFreePdofList.resize( mMyPdofHosts.size() );

        for ( moris::uint Ik=0; Ik < mMyPdofHosts.size(); Ik++ )
        {
            mFreePdofList( Ik ).resize( mMyPdofHosts( Ik )( 0 )->get_pdof_hosts_pdof_list().size() );
        }

        for ( moris::uint Ia=0; Ia < mMyPdofHosts.size(); Ia++ )
        {
            // Loop over all pdof hosts and get their number of (free) pdofs
            for ( moris::uint Ik=0; Ik < mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list().size(); Ik++ )
            {
                uint tNumPdofs = 0;

                for ( moris::uint Ii=0; Ii < mMyPdofHosts( Ia ).size(); Ii++ )
                {
                    tNumPdofs = tNumPdofs + mMyPdofHosts( Ia )( Ii )->get_pdof_hosts_pdof_list()( Ik ).size();
                }

                mFreePdofList( Ia )( Ik ).reserve( tNumPdofs );
            }
        }

        for ( moris::uint Ia=0; Ia < mMyPdofHosts.size(); Ia++ )
        {
            // Loop over all pdof hosts and get their number of (free) pdofs
            for ( moris::uint Ik=0; Ik < mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list().size(); Ik++ )
            {
                for ( moris::uint Ii=0; Ii < mMyPdofHosts( Ia ).size(); Ii++ )
                {
                    mFreePdofList( Ia )( Ik ).append( mMyPdofHosts( Ia )( Ii )->get_pdof_hosts_pdof_list()( Ik ) );
                }
            }
        }
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::create_my_list_of_adof_ids()
    {
        {
            // Get MAX number of pdofs for this equation object
            moris::uint tNumMyPdofs = mFreePdofs.size();

            // Loop over all pdofs to count their adofs
            moris::uint tNumMyAdofs = 0;
            for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
            {
                // Get Number of adofs cooresponding to this pdof
                moris::uint tNumAdofForThisPdof = ( mFreePdofs( Ij )->mAdofIds ).numel();
                tNumMyAdofs = tNumMyAdofs + tNumAdofForThisPdof;
            }

            // Temporary matrix for adofs Ids
            Matrix< DDSMat > tNonUniqueAdofIds( tNumMyAdofs, 1 );

            moris::uint tAdofPosCounter = 0;

            // Loop over all pdofs to get their adofs and put them into a unique list
            for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
            {
                tNonUniqueAdofIds ( {tAdofPosCounter, tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).numel() -1 }, { 0, 0} ) = mFreePdofs( Ij )->mAdofIds.matrix_data();

                // Add number if these adofs to number of assembled adofs
                tAdofPosCounter = tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).numel();
            }

            // make list of unique Ids
            moris::unique( tNonUniqueAdofIds, mUniqueAdofList );
        }

        //---------------------------------------------------------------------------

        mUniqueAdofTypeList.resize( mFreePdofList.size() );

        for ( moris::uint Ik=0; Ik < mFreePdofList.size(); Ik++ )
        {
            mUniqueAdofTypeList( Ik ).resize( mFreePdofList( Ik ).size() );
        }

        for ( moris::uint Ia=0; Ia < mFreePdofList.size(); Ia++ )
        {
            for ( moris::uint Ik=0; Ik < mFreePdofList( Ia ).size(); Ik++ )
            {
                // Get MAX number of pdofs for this equation object
                moris::uint tNumMyPdofs = mFreePdofList( Ia )( Ik ).size();

                // Loop over all pdofs to count their adofs
                moris::uint tNumMyAdofs = 0;
                for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
                {
                    // Get Number of adofs cooresponding to this pdof
                    moris::uint tNumAdofForThisPdof = ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel();
                    tNumMyAdofs += tNumAdofForThisPdof;
                }

                // Temporary matrix for adofs Ids
                Matrix< DDSMat > tNonUniqueAdofIds( tNumMyAdofs, 1 );

                moris::uint tAdofPosCounter = 0;

                // Loop over all pdofs to get their adofs and put them into a unique list
                for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
                {
                    tNonUniqueAdofIds ( { tAdofPosCounter, tAdofPosCounter + ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel() -1 }, { 0, 0 } )
                            = mFreePdofList( Ia )( Ik )( Ij )->mAdofIds.matrix_data();

                    // Add number if these adofs to number of assembled adofs
                    tAdofPosCounter += ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel();
                }

                // make list of unique Ids
                moris::unique( tNonUniqueAdofIds, mUniqueAdofTypeList( Ia )( Ik ) );
            }
        }
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::set_unique_adof_map()
    {
        // Loop over all unique adofs of this equation object
        for ( moris::uint Ii = 0; Ii < mUniqueAdofList.numel(); Ii++ )
        {
            mUniqueAdofMap[ mUniqueAdofList( Ii, 0 ) ] = Ii;
        }

        //----------------------------------------------------------

        moris::uint tNumAdofLists = mUniqueAdofTypeList.size();

        mUniqueAdofMapList.resize( tNumAdofLists );

        for ( moris::uint Ij = 0; Ij < tNumAdofLists; Ij++ )
        {
            //Get number of unique adofs of this equation object
            moris::uint tNumUniqueAdofsTypes = mUniqueAdofTypeList( Ij ).size();

            mUniqueAdofMapList( Ij ).resize( tNumUniqueAdofsTypes );

            // Loop over all unique adofs tpes of this equation object
            for ( moris::uint Ii = 0; Ii < tNumUniqueAdofsTypes; Ii++ )
            {
                moris::uint tNumUniqueAdofs = mUniqueAdofTypeList( Ij )( Ii ).numel();

                for ( moris::uint Ik = 0; Ik < tNumUniqueAdofs; Ik++ )
                {
                    mUniqueAdofMapList( Ij )( Ii )[ mUniqueAdofTypeList( Ij )( Ii )( Ik, 0 ) ] = Ik;
                }
            }
        }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::build_PADofMap( Matrix< DDRMat > & aPADofMap )
    {
         //Get number of unique adofs of this equation object
         moris::uint tNumUniqueAdofs = mUniqueAdofList.numel();

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
             for ( moris::uint Ik = 0; Ik < tPdof->mAdofIds.numel(); Ik++ )
             {
                 // Getting tPADofMap column entry for the corresponding value
                 moris::uint tColumnPos = mUniqueAdofMap[ tPdof->mAdofIds( Ik, 0 ) ];

                 // Insert value into pdof-adof-map
                 aPADofMap( Ii, tColumnPos ) = ( mFreePdofs( Ii )->mTmatrix)( Ik, 0 );
             }
         }
     }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::build_PADofMap_list( Cell< Cell< Matrix< DDRMat > > > & aPADofMap )
    {
        aPADofMap.resize( mUniqueAdofTypeList.size() );

        for ( moris::uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
        {
            moris::uint tNumUniqueAdofsTypes = mUniqueAdofTypeList( Ik ).size();

            aPADofMap( Ik ).resize( tNumUniqueAdofsTypes, Matrix<DDRMat> ( 0, 0 ) );
        }

        for ( moris::uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
        {
            // Loop over all adof types of this equation object
            for ( moris::uint Ij = 0; Ij < mUniqueAdofTypeList( Ik ).size(); Ij++ )
            {
                //Get number of unique adofs of this equation object
                moris::uint tNumUniqueAdofs = mUniqueAdofTypeList( Ik )( Ij ).numel();

                MORIS_ASSERT( tNumUniqueAdofs != 0,"Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. ");

                // Get MAX number of pdofs for this equation object
                moris::uint tNumMyPdofs = mFreePdofList( Ik )( Ij ).size();

                MORIS_ASSERT( tNumMyPdofs != 0,"Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. ");

                aPADofMap( Ik )( Ij ).set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

                // Loop over all pdofs of this equation object
                for ( moris::uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
                {
                    auto tPdof = mFreePdofList( Ik )( Ij )( Ii );

                    // Loop over all adof Ids of this pdof
                    for ( moris::uint Ib = 0; Ib < tPdof->mAdofIds.numel(); Ib++ )
                    {
                        // Getting tPADofMap column entry for the corresponding value
                        moris::uint tColumnPos = mUniqueAdofMapList( Ik )( Ij )[ tPdof->mAdofIds( Ib, 0 ) ];

                        // Insert value into pdof-adof-map
                        aPADofMap( Ik )( Ij )( Ii, tColumnPos ) = ( mFreePdofList( Ik )( Ij )( Ii )->mTmatrix)( Ib, 0 );
                    }
                }
            }
        }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::build_PADofMap_1( Matrix< DDRMat > & aPADofMap )
    {
        Cell< Cell< Matrix< DDRMat > > > tPADofMapList;

        // get list of all T-Matrices
        this->build_PADofMap_list( tPADofMapList );

        //get list of requested dof tpes
        moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes =  mEquationBlock->get_model_solver_interface()
                                                                              ->get_solver_interface()
                                                                              ->get_requested_dof_types();

        // initialize column and row counter
        uint tNumColCounter = 0;
        uint tNumRowCounter = 0;

        for ( moris::uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
        {
            // Loop over all requested dof types and get total number of cols and rows
            for ( moris::uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get index corresponding to this dof type
                moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()
                                                          ->get_dof_manager()
                                                          ->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                MORIS_ASSERT( tDofTypeIndex != -1,"build_PADofMap_1(), Index for this dof type does not exist in map");

                tNumColCounter += tPADofMapList( Ii )( tDofTypeIndex ).n_cols();
                tNumRowCounter += tPADofMapList( Ii )( tDofTypeIndex ).n_rows();
            }
        }

        aPADofMap.set_size( tNumRowCounter, tNumColCounter, 0.0 );

        // re-initialize column and row counter
        tNumColCounter = 0;
        tNumRowCounter = 0;

        for ( moris::uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
        {
            // Loop over all requested dof types and insert T-matrices into requested T-matrix
            for ( moris::uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get index corresponding to this dof type
                moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()
                                                          ->get_dof_manager()
                                                          ->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                uint tNumCols = tPADofMapList( Ii )( tDofTypeIndex ).n_cols();
                uint tNumRows = tPADofMapList( Ii )( tDofTypeIndex ).n_rows();

                aPADofMap( { tNumRowCounter, tNumRowCounter + tNumRows -1 }, { tNumColCounter, tNumColCounter+ tNumCols - 1 } ) += tPADofMapList( Ii )( tDofTypeIndex ).matrix_data();

                tNumColCounter += tNumCols;
                tNumRowCounter += tNumRows;
            }
        }
    }

//-------------------------------------------------------------------------------------------------
    moris_index Equation_Object::get_node_index( const moris_index aElementLocalNodeIndex ) const
    {
        return mNodeObj( 0 )( aElementLocalNodeIndex )->get_index();
    }

//-------------------------------------------------------------------------------------------------
    void Equation_Object::get_equation_obj_dof_ids( Matrix< DDSMat > & aEqnObjAdofId )
    {
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  mEquationBlock->get_model_solver_interface()
                                                                         ->get_solver_interface()
                                                                         ->get_requested_dof_types();

        Dof_Manager * tDofManager =  mEquationBlock->get_model_solver_interface()->get_dof_manager();

        uint tCounter = 0;
        for ( moris::uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
        {
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                moris::sint tDofTypeIndex = tDofManager->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                tCounter += mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();
            }
        }

        aEqnObjAdofId.set_size( tCounter, 1, -1 );

        tCounter = 0;
        for ( moris::uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
        {
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                moris::sint tDofTypeIndex = tDofManager->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                aEqnObjAdofId( { tCounter, tCounter + mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel() - 1 }, { 0, 0 }) = mUniqueAdofTypeList( Ii )( tDofTypeIndex ).matrix_data();

                tCounter += mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();
            }
        }

        MORIS_ASSERT( aEqnObjAdofId.min() != -1, "Equation_Onject::get_equation_obj_dof_ids(), Error while returning adof ids for type" );

//        print (aEqnObjAdofId,"aEqnObjAdofId");
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::get_egn_obj_jacobian( Matrix< DDRMat > & aEqnObjMatrix,
                                                Dist_Vector      * aSolutionVector )
    {
        mSolVec = aSolutionVector;

        // compute jacobin
        this->compute_jacobian();


        // build T-matrix
        Matrix< DDRMat > tTMatrix;
        this->build_PADofMap_1( tTMatrix );

//        print( tTMatrix,"tTMatrix");
//        print( mEquationBlock->get_jacobian(),"tJacobian");



        // project pdof resdiual to adof residual
        aEqnObjMatrix = trans( tTMatrix ) * mEquationBlock->get_jacobian() * tTMatrix;

//        print(aEqnObjMatrix,"aEqnObjMatrix");

    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::get_equation_obj_residual( Matrix< DDRMat > & aEqnObjRHS, Dist_Vector * aSolutionVector )
    {
        mSolVec = aSolutionVector;

        this->compute_residual();

        Matrix< DDRMat > tElementalResidual = mEquationBlock->get_residual();

//        print( mEquationBlock->get_residual(),"Residual");

        this->add_staggered_contribution_to_residual( tElementalResidual );

        Matrix< DDRMat > tTMatrix;
        this->build_PADofMap_1( tTMatrix );

        aEqnObjRHS = trans( tTMatrix ) * tElementalResidual;

//        print(aEqnObjRHS, "aEqnObjRHS");

    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::add_staggered_contribution_to_residual( Matrix< DDRMat > & aElementResidual )
    {
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tSecDofTypes =  mEquationBlock->get_secundary_dof_types();

       for( auto tSecDofTypes : tSecDofTypes )
       {
            moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes =  mEquationBlock->get_requested_dof_types();

            // combined master slave index
            sint tSecDofIndex = mEquationBlock->get_dof_index_for_type( tSecDofTypes( 0 ), mtk::Master_Slave::MASTER );

            if( tSecDofIndex != -1 )
            {
                for( auto tDofTypes : tRequestedDofTypes )
                {
                    sint tDofIndex = mEquationBlock->get_dof_index_for_type( tDofTypes, mtk::Master_Slave::MASTER );

                    if( tDofIndex != -1 )
                    {
                        uint tStartRow = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 0 );
                        uint tEndRow   = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 1 );

                        uint tStartCol = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                        uint tEndCol   = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                        // seperate master slave index
                        sint tSeperateSecDofIndex = mEquationBlock->get_dof_index_for_type_1( tSecDofTypes( 0 ), mtk::Master_Slave::MASTER );

                        // get the ith dof type group
                        moris::Cell< MSI::Dof_Type > tDofTypeGroup = mEquationBlock->mMasterDofTypes( tSeperateSecDofIndex );

                        // get the pdof values for the ith dof type group
                        Cell< Matrix< DDRMat > > tCoeff_Original;
                        Matrix< DDRMat > tCoeff;

                        this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::MASTER );

                        // reshape tCoeffs into the order the cluster expects them
                        this->reshape_pdof_values( tCoeff_Original, tCoeff );

                        Matrix< DDRMat > tCoeff1( tCoeff.numel(),1, 0.0);
                        //fixme get rid of for loop
                        uint tCounter=0;
                        for( uint Ik = 0; Ik<tCoeff.n_cols(); Ik++)
                        {
                            for( uint Ii = 0; Ii<tCoeff.n_rows(); Ii++)
                            {
                                tCoeff1(tCounter++) = tCoeff(Ii,Ik);
                            }
                        }

//                        print(aElementResidual,"aElementResidual");
//                        print(tCoeff1,"tCoeff1");
//                        print(mEquationBlock->get_jacobian(),"mEquationBlock->get_jacobian()");
                        aElementResidual( { tStartRow, tEndRow },{ 0, 0 } )
                               -= mEquationBlock->get_jacobian()( { tStartRow, tEndRow }, { tStartCol, tEndCol } ) * tCoeff1;
                    }

                    tDofIndex = mEquationBlock->get_dof_index_for_type( tDofTypes, mtk::Master_Slave::SLAVE );

                    if( tDofIndex != -1 )
                    {
                        uint tStartRow = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 0 );
                        uint tEndRow   = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 1 );

                        uint tStartCol = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                        uint tEndCol   = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                        // seperate master slave index
                        sint tSeperateSecDofIndex = mEquationBlock->get_dof_index_for_type_1( tSecDofTypes( 0 ), mtk::Master_Slave::SLAVE );

                        // get the ith dof type group
                        moris::Cell< MSI::Dof_Type > tDofTypeGroup = mEquationBlock->mSlaveDofTypes( tSeperateSecDofIndex );

                        // get the pdof values for the ith dof type group
                        Cell< Matrix< DDRMat > > tCoeff_Original;
                        Matrix< DDRMat > tCoeff;

                        this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

                        // reshape tCoeffs into the order the cluster expects them
                        this->reshape_pdof_values( tCoeff_Original, tCoeff );

                        aElementResidual( { tStartRow, tEndRow },{ 0, 0 } )
                               -= mEquationBlock->get_jacobian()( { tStartRow, tEndRow }, { tStartCol, tEndCol } ) * tCoeff;
                    }
                }
            }

            // combined master slave index
            tSecDofIndex = mEquationBlock->get_dof_index_for_type( tSecDofTypes( 0 ), mtk::Master_Slave::SLAVE );

            if( tSecDofIndex != -1 )
            {
                for( auto tDofTypes : tRequestedDofTypes )
                {
                    sint tDofIndex = mEquationBlock->get_dof_index_for_type( tDofTypes, mtk::Master_Slave::MASTER );

                    if( tDofIndex != -1 )
                    {
                        uint tStartRow = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 0 );
                        uint tEndRow   = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 1 );

                        uint tStartCol = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                        uint tEndCol   = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                        // seperate master slave index
                        sint tSeperateSecDofIndex = mEquationBlock->get_dof_index_for_type_1( tSecDofTypes( 0 ), mtk::Master_Slave::MASTER );

                        // get the ith dof type group
                        moris::Cell< MSI::Dof_Type > tDofTypeGroup = mEquationBlock->mMasterDofTypes( tSeperateSecDofIndex );

                        // get the pdof values for the ith dof type group
                        Cell< Matrix< DDRMat > > tCoeff_Original;
                        Matrix< DDRMat > tCoeff;

                        this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::MASTER );

                        // reshape tCoeffs into the order the cluster expects them
                        this->reshape_pdof_values( tCoeff_Original, tCoeff );

                        aElementResidual( { tStartRow, tEndRow },{ 0, 0 } )
                               -= mEquationBlock->get_jacobian()( { tStartRow, tEndRow }, { tStartCol, tEndCol } ) * tCoeff;
                    }

                    tDofIndex = mEquationBlock->get_dof_index_for_type( tDofTypes, mtk::Master_Slave::SLAVE );

                    if( tDofIndex != -1 )
                    {
                        uint tStartRow = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 0 );
                        uint tEndRow   = mEquationBlock->mResDofAssemblyMap( tDofIndex )( 0, 1 );

                        uint tStartCol = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                        uint tEndCol   = mEquationBlock->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                        // seperate master slave index
                        sint tSeperateSecDofIndex = mEquationBlock->get_dof_index_for_type_1( tSecDofTypes( 0 ), mtk::Master_Slave::SLAVE );

                        // get the ith dof type group
                        moris::Cell< MSI::Dof_Type > tDofTypeGroup = mEquationBlock->mSlaveDofTypes( tSeperateSecDofIndex );

                        // get the pdof values for the ith dof type group
                        Cell< Matrix< DDRMat > > tCoeff_Original;
                        Matrix< DDRMat > tCoeff;

                        this->get_my_pdof_values( tDofTypeGroup, tCoeff_Original, mtk::Master_Slave::SLAVE );

                        // reshape tCoeffs into the order the cluster expects them
                        this->reshape_pdof_values( tCoeff_Original, tCoeff );

                        aElementResidual( { tStartRow, tEndRow },{ 0, 0 } )
                               -= mEquationBlock->get_jacobian()( { tStartRow, tEndRow }, { tStartCol, tEndCol } ) * tCoeff;
                    }
                }
            }
       }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::compute_my_pdof_values( )
    {
        Matrix< DDRMat > tTMatrix;

        // build T-matrix
        this->build_PADofMap( tTMatrix );

        Matrix< DDRMat > tMyValues;

        // Extract this equation objects adof values from solution vector
        mSolVec->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

        // multiply t_matrix with adof values to get pdof values
        mPdofValues = tTMatrix * tMyValues;

        this->set_vector_entry_number_of_pdof();             // FIXME should not be in MSI. Should be in FEM
    }

//-------------------------------------------------------------------------------------------------

//    void Equation_Object::get_my_pdof_values( const moris::Cell< enum Dof_Type > & aRequestedDofTypes,
//                                                    Cell< Matrix< DDRMat > >     & aRequestedPdofValues )
//    {
//        // Initialize list which contains the maximal number of time levels per dof type
//        Matrix< DDSMat > tTimeLevelsPerDofType( aRequestedDofTypes.size(), 1, -1 );
//
//        moris::sint tCounter = 0;
//
//        // Loop over requested dof types
//        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
//        {
//            // Loop over all elemental pdof hosts
//            for ( moris::uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
//            {
//                // Get dof type index
//                moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()->get_dof_manager()
//                                                                 ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );
//
//                MORIS_ERROR( mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) !=0,
//                        "Equation_Object::get_my_pdof_values: talk with Mathias about this");                         //FIXME delete this error after a closer look
//
//                // get number of time levels for this dof type
//                moris::sint tNumTimeLevels = mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex );
//                tCounter = tCounter + tNumTimeLevels;
//
//                // Add maximal value of time levels to list
//                tTimeLevelsPerDofType( Ii, 0 ) = std::max( tTimeLevelsPerDofType( Ii, 0 ), tNumTimeLevels );
//            }
//            MORIS_ASSERT( tTimeLevelsPerDofType( Ii, 0 ) > -1, "Equation_Object::get_my_pdof_values: no time levels exist on this dof type on element %-5i", mEqnObjInd );
//        }
//        // Set size matrix for requested pdof values
//        aRequestedPdofValues.resize( tCounter, 1 );
//
//        moris::sint tCounter_2 = 0;
//
//        // Loop over requested dof types
//        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
//        {
//            // Get maximal Number of time levels on this pdof type
//            moris::sint tMaxTimeLevelsOnDofType = tTimeLevelsPerDofType( Ii, 0 );
//
//            // Loop over this pdofs time levels
//            for ( moris::sint Ia = 0; Ia < tMaxTimeLevelsOnDofType; Ia++ )
//            {
//                // Loop over all elemental pdof hosts
//                for ( moris::uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
//                {
//                    // Get dof type index
//                    moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()->get_dof_manager()
//                                                                     ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );
//
//                    // Check if number if time levels on this dof type is smaller than maximal number of time levels on dof type
//                    if ( (sint)mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) == tMaxTimeLevelsOnDofType )
//                    {
//                        // get pointer list all time pdofs on this pdof type
//                        moris::Cell< Pdof* > tPdofTimeList = mMyPdofHosts( Ik )->get_pdof_time_list( tDofTypeIndex );
//
//                        // get entry number of this pdof in the elemental pdof value vector
//                        moris::uint mElementalSolVecEntry = tPdofTimeList( Ia )->mElementalSolVecEntry;
//
//                        // Put this pdof value into the requested pdof vector
//                        aRequestedPdofValues( tCounter_2++, 0 ) = mPdofValues( mElementalSolVecEntry , 0 );
//                    }
//                }
//            }
//        }
//    }

    void Equation_Object::get_my_pdof_values( const moris::Cell< enum Dof_Type > & aRequestedDofTypes,
                                                    Cell< Matrix< DDRMat > >     & aRequestedPdofValues,
                                                    mtk::Master_Slave              aIsMaster )
    {
        uint tIsMaster = 0;

        switch ( aIsMaster )
        {
            case ( mtk::Master_Slave::MASTER ):
            {
                 tIsMaster = 0;
                 break;
            }
            case( mtk::Master_Slave::SLAVE ):
            {
                tIsMaster = 1;
                break;
            }
            default:
            {
                MORIS_ERROR(false, "Equation_Object::get_my_pdof_values - can only be MASTER or SLAVE");
            }
        }

        // Initialize list which contains the maximal number of time levels per dof type
        Matrix< DDSMat > tTimeLevelsPerDofType( aRequestedDofTypes.size(), 1, -1 );

        aRequestedPdofValues.resize( aRequestedDofTypes.size() );

        moris::sint tCounter = 0;

        // Loop over requested dof types
        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        {
            tCounter = 0;

            // Loop over all elemental pdof hosts
            for ( moris::uint Ik = 0; Ik < mMyPdofHosts( tIsMaster ).size(); Ik++ )
            {
                // Get dof type index

                moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()->get_dof_manager()
                                                          ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                MORIS_ASSERT( mMyPdofHosts( tIsMaster )( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) !=0,
                        "Equation_Object::get_my_pdof_values: talk with Mathias about this");                         //FIXME delete this error after a closer look

                // get number of time levels for this dof type
                moris::sint tNumTimeLevels = mMyPdofHosts( tIsMaster )( Ik )->get_num_time_levels_of_type( tDofTypeIndex );
                tCounter = tCounter + tNumTimeLevels;

                // Add maximal value of time levels to list
                tTimeLevelsPerDofType( Ii, 0 ) = std::max( tTimeLevelsPerDofType( Ii, 0 ), tNumTimeLevels );
            }
            MORIS_ASSERT( tTimeLevelsPerDofType( Ii, 0 ) > -1, "Equation_Object::get_my_pdof_values: no time levels exist on this dof type on element %-5i", mEqnObjInd );

            // Set size matrix for requested pdof values
            aRequestedPdofValues( Ii ).resize( tCounter, 1 );
        }

        moris::sint tCounter_2 = 0;

        // Loop over requested dof types
        for ( moris::uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        {
            tCounter_2 = 0;
            // Get maximal Number of time levels on this pdof type
            moris::sint tMaxTimeLevelsOnDofType = tTimeLevelsPerDofType( Ii, 0 );

            // Loop over this pdofs time levels
            for ( moris::sint Ia = 0; Ia < tMaxTimeLevelsOnDofType; Ia++ )
            {
                // Loop over all elemental pdof hosts
                for ( moris::uint Ik = 0; Ik < mMyPdofHosts( tIsMaster ).size(); Ik++ )
                {
                    // Get dof type index

                    moris::sint tDofTypeIndex = mEquationBlock->get_model_solver_interface()->get_dof_manager()
                                                              ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                    // Check if number if time levels on this dof type is smaller than maximal number of time levels on dof type
                    if ( (sint)mMyPdofHosts( tIsMaster )( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) == tMaxTimeLevelsOnDofType )
                    {
                        // get pointer list all time pdofs on this pdof type
                        moris::Cell< Pdof* > tPdofTimeList = mMyPdofHosts( tIsMaster )( Ik )->get_pdof_time_list( tDofTypeIndex );

                        // get entry number of this pdof in the elemental pdof value vector
                        moris::uint tElementalSolVecEntry = tPdofTimeList( Ia )->mElementalSolVecEntry;

                        // Put this pdof value into the requested pdof vector
                        aRequestedPdofValues( Ii )( tCounter_2++, 0 ) = mPdofValues( tElementalSolVecEntry , 0 );
                    }
                }
            }
        }
    }

//-------------------------------------------------------------------------------------------------

    void Equation_Object::reshape_pdof_values( const Cell< Matrix< DDRMat > > & aPdofValues,
                                                     Matrix< DDRMat >         & aReshapedPdofValues )
    {
        MORIS_ASSERT( aPdofValues.size() != 0, "Cluster::reshape_pdof_values(), pdof value vector is empty");

        uint tCols = aPdofValues.size();
        uint tRows = aPdofValues( 0 ).numel();

        aReshapedPdofValues.set_size( tRows, tCols );

        for( uint Ik = 0; Ik < tCols; Ik++ )
        {
            aReshapedPdofValues( { 0, tRows - 1 }, { Ik, Ik } ) = aPdofValues( Ik ).matrix_data();
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
