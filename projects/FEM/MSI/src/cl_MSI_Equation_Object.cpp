/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Object.cpp
 *
 */

#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "fn_isfinite.hpp"

namespace moris
{
    namespace MSI
    {

        //-------------------------------------------------------------------------------------------------

        uint
        Equation_Object::get_max_pdof_hosts_ind()
        {
            auto tMaxPdofHostsInd = mNodeObj( 0 )( 0 )->get_index();

            // Loop over all node obj. get the maximal node index.
            for ( uint Ik = 0; Ik < mNodeObj.size(); Ik++ )
            {
                for ( uint Ii = 0; Ii < mNodeObj( Ik ).size(); Ii++ )
                {
                    tMaxPdofHostsInd = std::max( tMaxPdofHostsInd, mNodeObj( Ik )( Ii )->get_index() );
                }
            }
            return (uint)tMaxPdofHostsInd;
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::create_my_pdof_hosts(
                const uint                 aNumUsedDofTypes,
                const Matrix< DDSMat >&    aPdofTypeMap,
                const Matrix< DDUMat >&    aTimePerDofType,
                Vector< Pdof_Host* >& aPdofHostList )
        {
            // Resize list containing this equations objects pdof hosts set
            mNumPdofSystems = mNodeObj.size();
            mMyPdofHosts.resize( mNumPdofSystems );    // FIXME: Add ghost and element numbers

            // Loop over pdof systems. Is one except for double sided clusters
            for ( uint Ik = 0; Ik < mNumPdofSystems; Ik++ )
            {
                uint tNumMyPdofHosts = mNodeObj( Ik ).size();

                // Resize list containing this equations objects pdof hosts
                mMyPdofHosts( Ik ).resize( tNumMyPdofHosts, nullptr );

                // Loop over all nodes of this element, creating new pdof hosts if not existing yet.
                for ( uint Ii = 0; Ii < mNodeObj( Ik ).size(); Ii++ )
                {
                    // Save node id of node Ii in temporary variable for more clarity.
                    // auto tNodeID = mNodeObj( Ik )( Ii )->get_id();
                    auto tNodeID = mNodeObj( Ik )( Ii )->get_index();

                    // check if pdof host corresponding to this node exists.
                    if ( aPdofHostList( tNodeID ) == NULL )
                    {
                        // If node does not exist, create new pdof host.
                        aPdofHostList( tNodeID ) = new Pdof_Host( aNumUsedDofTypes, mNodeObj( Ik )( Ii ) );
                    }

                    // Add pointer to pdof host to the list containing this equation objects pdof hosts.
                    mMyPdofHosts( Ik )( Ii ) = aPdofHostList( tNodeID );

                    uint tNumUsedDofTypes = mEquationSet->get_unique_leader_follower_dof_type_list()( Ik ).size();

                    // FIXME: rewrite this function
                    for ( uint iDofType = 0; iDofType < tNumUsedDofTypes; iDofType++ )
                    {
                        mMyPdofHosts( Ik )( Ii )->set_pdof_type(
                                mEquationSet->get_unique_leader_follower_dof_type_list()( Ik )( iDofType ),
                                aTimePerDofType,
                                aNumUsedDofTypes,
                                aPdofTypeMap );
                    }
                }
            }

            // FIXME: add element
            // FIXME: return pointer to pdofs
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::create_my_pdof_list()
        {
            // init number of free pdof counter
            uint tNumMyFreePdofs = 0;

            // loop over all pdof hosts and get their number of free pdofs
            for ( uint Ik = 0; Ik < mNumPdofSystems; Ik++ )
            {
                // get number of pdof hosts corresponding to this equation object
                uint tNumMyPdofHosts = mMyPdofHosts( Ik ).size();
                for ( uint Ii = 0; Ii < tNumMyPdofHosts; Ii++ )
                {
                    // add number of free pdofs to counter
                    tNumMyFreePdofs = tNumMyFreePdofs + mMyPdofHosts( Ik )( Ii )->get_num_pdofs();
                }
            }

            // set size of vector containing this equation objects free pdofs
            mFreePdofs.reserve( tNumMyFreePdofs );

            // loop over pdof systems. Is one except for double sided clusters
            for ( uint Ia = 0; Ia < mNumPdofSystems; Ia++ )
            {
                uint tNumMyPdofHosts = mMyPdofHosts( Ia ).size();

                // loop over all pdof types. Ask the first pdof host for the number of pdof types
                uint tNumPdofTypes = ( mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list() ).size();
                for ( uint Ij = 0; Ij < tNumPdofTypes; Ij++ )
                {
                    // loop over all time levels for this dof type
                    uint tNumTimeLevels = mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list()( Ij ).size();
                    for ( uint Ii = 0; Ii < tNumTimeLevels; Ii++ )
                    {
                        // loop over all pdof hosts and dof types. Appending the pdof pointers to the pdof list of this equation object
                        for ( uint Ik = 0; Ik < tNumMyPdofHosts; Ik++ )
                        {
                            // append all time levels of this pdof type
                            mFreePdofs.push_back( ( mMyPdofHosts( Ia )( Ik )->get_pdof_hosts_pdof_list() )(Ij)( Ii ) );
                        }
                    }
                }
            }

            //----------------------------------------------------------------------------------------------------------

            // Ask the first pdof host for the number of pdof types //FIXME
            mFreePdofList.resize( mNumPdofSystems );

            for ( uint Ia = 0; Ia < mNumPdofSystems; Ia++ )
            {
                mFreePdofList( Ia ).resize( mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list().size() );

                // Loop over all pdof hosts and get their number of (free) pdofs
                uint tNumPdofHosts = mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list().size();
                for ( uint Ik = 0; Ik < tNumPdofHosts; Ik++ )
                {
                    uint tNumPdofs    = 0;
                    uint tNumFreePdof = mMyPdofHosts( Ia ).size();
                    for ( uint Ii = 0; Ii < tNumFreePdof; Ii++ )
                    {
                        tNumPdofs = tNumPdofs + mMyPdofHosts( Ia )( Ii )->get_pdof_hosts_pdof_list()( Ik ).size();
                    }
                    mFreePdofList( Ia )( Ik ).reserve( tNumPdofs );
                }
            }

            for ( uint Ia = 0; Ia < mNumPdofSystems; Ia++ )
            {
                uint tNumMyPdofHosts = mMyPdofHosts( Ia ).size();

                // Loop over all pdof hosts and get their number of (free) pdofs
                uint tNumPdofHosts = mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list().size();
                for ( uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
                {
                    // Loop over all time levels for this dof type
                    uint tNumTimeLevels = mMyPdofHosts( Ia )( 0 )->get_pdof_hosts_pdof_list()( Ij ).size();
                    for ( uint Ii = 0; Ii < tNumTimeLevels; Ii++ )
                    {
                        for ( uint Ik = 0; Ik < tNumMyPdofHosts; Ik++ )
                        {
                            // mFreePdofList( Ia )( Ik ).append( mMyPdofHosts( Ia )( Ii )->get_pdof_hosts_pdof_list()( Ik ) );
                            // Append all time levels of this pdof type
                            mFreePdofList( Ia )( Ij ).push_back( ( mMyPdofHosts( Ia )( Ik )->get_pdof_hosts_pdof_list() )(Ij)( Ii ) );
                        }
                    }
                }
            }

            // free pdof type list created, set flag to true
            mFreePdofListFlag = true;
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::create_my_list_of_adof_ids()
        {
            {
                // Get MAX number of pdofs for this equation object
                uint tNumMyPdofs = mFreePdofs.size();

                // Loop over all pdofs to count their adofs
                uint tNumMyAdofs = 0;
                for ( uint Ij = 0; Ij < tNumMyPdofs; Ij++ )
                {
                    // Get Number of adofs corresponding to this pdof
                    uint tNumAdofForThisPdof = ( mFreePdofs( Ij )->mAdofIds ).numel();
                    tNumMyAdofs              = tNumMyAdofs + tNumAdofForThisPdof;
                }

                // Temporary matrix for adofs Ids
                Matrix< DDSMat > tNonUniqueAdofIds( tNumMyAdofs, 1 );

                uint tAdofPosCounter = 0;

                // Loop over all pdofs to get their adofs and put them into a unique list
                for ( uint Ij = 0; Ij < tNumMyPdofs; Ij++ )
                {
                    if ( ( mFreePdofs( Ij )->mAdofIds ).numel() == 0 ) continue; 
			tNonUniqueAdofIds( { tAdofPosCounter, tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).numel() - 1 }, { 0, 0 } ) =
                            mFreePdofs( Ij )->mAdofIds.matrix_data();

                    // Add number if these adofs to number of assembled adofs
                    tAdofPosCounter = tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).numel();
                }

                // make list of unique Ids
                moris::unique( tNonUniqueAdofIds, mUniqueAdofList );
            }

            //---------------------------------------------------------------------------

            mUniqueAdofTypeList.resize( mFreePdofList.size() );

            for ( uint Ik = 0; Ik < mFreePdofList.size(); Ik++ )
            {
                mUniqueAdofTypeList( Ik ).resize( mFreePdofList( Ik ).size() );
            }

            for ( uint Ia = 0; Ia < mFreePdofList.size(); Ia++ )
            {
                for ( uint Ik = 0; Ik < mFreePdofList( Ia ).size(); Ik++ )
                {
                    // Get MAX number of pdofs for this equation object
                    uint tNumMyPdofs = mFreePdofList( Ia )( Ik ).size();

                    // Loop over all pdofs to count their adofs
                    uint tNumMyAdofs = 0;
                    for ( uint Ij = 0; Ij < tNumMyPdofs; Ij++ )
                    {
                        // Get Number of adofs corresponding to this pdof
                        uint tNumAdofForThisPdof = ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel();
                        tNumMyAdofs += tNumAdofForThisPdof;
                    }

                    // Temporary matrix for adofs Ids
                    Matrix< DDSMat > tNonUniqueAdofIds( tNumMyAdofs, 1 );

                    uint tAdofPosCounter = 0;

                    // Loop over all pdofs to get their adofs and put them into a unique list
                    for ( uint Ij = 0; Ij < tNumMyPdofs; Ij++ )
                    {
                    if ( (mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel() == 0 ) continue; 
			    uint tEndAddress = tAdofPosCounter + ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel() - 1;

                        tNonUniqueAdofIds( { tAdofPosCounter, tEndAddress } ) =
                                mFreePdofList( Ia )( Ik )( Ij )->mAdofIds.matrix_data();

                        // Add number if these adofs to number of assembled adofs
                        tAdofPosCounter += ( mFreePdofList( Ia )( Ik )( Ij )->mAdofIds ).numel();
                    }

                    // make list of unique Ids
                    moris::unique( tNonUniqueAdofIds, mUniqueAdofTypeList( Ia )( Ik ) );
                }
            }

            // unique adof type list created, set flag to true
            mUniqueAdofTypeListFlag = true;
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::set_unique_adof_map()
        {
            // Loop over all unique adofs of this equation object
            for ( uint Ii = 0; Ii < mUniqueAdofList.numel(); Ii++ )
            {
                mUniqueAdofMap[ mUniqueAdofList( Ii, 0 ) ] = Ii;
            }

            //----------------------------------------------------------

            uint tNumAdofLists = mUniqueAdofTypeList.size();

            mUniqueAdofMapList.resize( tNumAdofLists );

            for ( uint Ij = 0; Ij < tNumAdofLists; Ij++ )
            {
                // Get number of unique adofs of this equation object
                uint tNumUniqueAdofsTypes = mUniqueAdofTypeList( Ij ).size();

                mUniqueAdofMapList( Ij ).resize( tNumUniqueAdofsTypes );

                // Loop over all unique adofs types of this equation object
                for ( uint Ii = 0; Ii < tNumUniqueAdofsTypes; Ii++ )
                {
                    uint tNumUniqueAdofs = mUniqueAdofTypeList( Ij )( Ii ).numel();

                    for ( uint Ik = 0; Ik < tNumUniqueAdofs; Ik++ )
                    {
                        mUniqueAdofMapList( Ij )( Ii )[ mUniqueAdofTypeList( Ij )( Ii )( Ik, 0 ) ] = Ik;
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_PADofMap( Matrix< DDRMat >& aPADofMap )
        {
            // Get number of unique adofs of this equation object
            uint tNumUniqueAdofs = mUniqueAdofList.numel();

            //MORIS_ASSERT( tNumUniqueAdofs != 0,
            //        "Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. " );

            // Get MAX number of pdofs for this equation object
            uint tNumMyPdofs = mFreePdofs.size();

            MORIS_ASSERT( tNumMyPdofs != 0,
                    "Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. " );

            aPADofMap.set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

            // Loop over all pdofs of this equation object
            for ( uint iPDOF = 0; iPDOF < tNumMyPdofs; iPDOF++ )
            {
                auto tPdof = mFreePdofs( iPDOF );

                // Loop over all adof Ids of this pdof
                for ( uint Ik = 0; Ik < tPdof->mAdofIds.numel(); Ik++ )
                {
                    // Getting tPADofMap column entry for the corresponding value
                    uint tColumnPos = mUniqueAdofMap[ tPdof->mAdofIds( Ik, 0 ) ];

                    // Insert value into pdof-adof-map
                    aPADofMap( iPDOF, tColumnPos ) = ( mFreePdofs( iPDOF )->mTmatrix )( Ik, 0 );
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_PADofMap_list( Vector< Vector< Matrix< DDRMat > > >& aPADofMap )
        {
            aPADofMap.resize( mUniqueAdofTypeList.size() );

            for ( uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
            {
                uint tNumUniqueAdofsTypes = mUniqueAdofTypeList( Ik ).size();

                aPADofMap( Ik ).resize( tNumUniqueAdofsTypes, Matrix< DDRMat >( 0, 0 ) );
            }

            for ( uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
            {
                // Loop over all adof types of this equation object
                for ( uint Ij = 0; Ij < mUniqueAdofTypeList( Ik ).size(); Ij++ )
                {
                    MORIS_ASSERT( mUniqueAdofTypeListFlag,
                            "Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. " );

                    // Get number of unique adofs of this equation object
                    uint tNumUniqueAdofs = mUniqueAdofTypeList( Ik )( Ij ).numel();

                    // MORIS_ASSERT( tNumUniqueAdofs != 0,
                    //         "Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. ");

                    MORIS_ASSERT( mFreePdofListFlag,
                            "Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. " );

                    // Get MAX number of pdofs for this equation object
                    uint tNumMyPdofs = mFreePdofList( Ik )( Ij ).size();

                    // MORIS_ASSERT( tNumMyPdofs != 0,
                    //        "Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. ");

                    aPADofMap( Ik )( Ij ).set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

                    // Loop over all pdofs of this equation object
                    for ( uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
                    {
                        auto tPdof = mFreePdofList( Ik )( Ij )( Ii );

                        // Loop over all adof Ids of this pdof
                        for ( uint Ib = 0; Ib < tPdof->mAdofIds.numel(); Ib++ )
                        {
                            // Getting tPADofMap column entry for the corresponding value
                            uint tColumnPos = mUniqueAdofMapList( Ik )( Ij )[ tPdof->mAdofIds( Ib, 0 ) ];

                            // Insert value into pdof-adof-map
                            aPADofMap( Ik )( Ij )( Ii, tColumnPos ) = ( mFreePdofList( Ik )( Ij )( Ii )->mTmatrix )( Ib, 0 );
                        }
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_PADofMap_1( Matrix< DDRMat >& aPADofMap )
        {
            Vector< Vector< Matrix< DDRMat > > > tPADofMapList;

            // get list of all T-Matrices
            this->build_PADofMap_list( tPADofMapList );

            Vector< enum MSI::Dof_Type > tRequestedDofTypes;

            // get list of requested dof types
            if ( !mEquationSet->mIsStaggered )
            {
                tRequestedDofTypes = mEquationSet->get_requested_dof_types();
            }
            else
            {
                tRequestedDofTypes = mEquationSet->get_secondary_dof_types();
            }

            // initialize column and row counter
            uint tNumColCounter = 0;
            uint tNumRowCounter = 0;

            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                // Loop over all requested dof types and get total number of cols and rows
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get index corresponding to this dof type
                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    tNumColCounter += tPADofMapList( Ii )( tDofTypeIndex ).n_cols();
                    tNumRowCounter += tPADofMapList( Ii )( tDofTypeIndex ).n_rows();
                }
            }
            aPADofMap.set_size( tNumRowCounter, tNumColCounter, 0.0 );

            // re-initialize column and row counter
            tNumColCounter = 0;
            tNumRowCounter = 0;

            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                // Loop over all requested dof types and insert T-matrices into requested T-matrix
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get index corresponding to this dof type
                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    // get number of rows and columns
                    uint tNumCols = tPADofMapList( Ii )( tDofTypeIndex ).n_cols();
                    uint tNumRows = tPADofMapList( Ii )( tDofTypeIndex ).n_rows();

                    // tPADofMapList( Ii )( tDofTypeIndex ) is not empty
                    if ( tNumCols * tNumRows > 0 )
                    {
                        aPADofMap(
                                { tNumRowCounter, tNumRowCounter + tNumRows - 1 },
                                { tNumColCounter, tNumColCounter + tNumCols - 1 } ) +=
                                tPADofMapList( Ii )( tDofTypeIndex ).matrix_data();

                        tNumColCounter += tNumCols;
                        tNumRowCounter += tNumRows;
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_dPADofMap_list( Vector< Vector< Vector< Vector< Vector< real > > > > >&        adPADofMap,
                                               Vector< Vector< Vector< Vector< Vector< moris_id > > > > >&    adPADofMapInds )
        {
            // This function is used to create d(PADoFMap)/d\hat{x}, i.e., the derivative of the PADofMap with respect to the nodal coordinates.
            // Why do we need it? - Basis extension leads to the T-Matrix also posessing sensitivities with respect to nodal locations
            // Due to the averaging weights used for basis extension to get the volume fraction.

            // the nesting is such - outermost is the leader/follower system (bulk has only 1)
            // second is the DoF types (UX, UY etc.)
            // third is the pdofs of the type above associated with the equation object
            // 4th is the adofs associated with each pdof - the col index from the build PADofMap list func is used to index this
            // the first input argument stores the x,y,z components of the derivative per adof, and second stores the IG node IDs w.r.t.
            // which there are sensitivities.
            // the DDRMats will be num PDVs

            // get spatial dimensions
            Matrix< DDRMat > tCoords;
            mNodeObj( 0 )( 0 )->get_vertex_coords( tCoords );
            uint tNumSpatialDims = tCoords.numel();

            // Remove any results from a previous call.
            adPADofMap.clear();
            adPADofMapInds.clear();

            adPADofMap.resize( mUniqueAdofTypeList.size() );
            adPADofMapInds.resize( mUniqueAdofTypeList.size() );

            for ( uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
            {
                uint tNumUniqueAdofsTypes = mUniqueAdofTypeList( Ik ).size();

                adPADofMap( Ik ).resize( tNumUniqueAdofsTypes );
                adPADofMapInds( Ik ).resize( tNumUniqueAdofsTypes );
            }

            for ( uint Ik = 0; Ik < mUniqueAdofTypeList.size(); Ik++ )
            {
                // Loop over all adof types of this equation object
                for ( uint Ij = 0; Ij < mUniqueAdofTypeList( Ik ).size(); Ij++ )
                {
                    MORIS_ASSERT( mUniqueAdofTypeListFlag,
                            "Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. " );

                    // Get number of unique adofs of this equation object
                    uint tNumUniqueAdofs = mUniqueAdofTypeList( Ik )( Ij ).numel();

                    // MORIS_ASSERT( tNumUniqueAdofs != 0,
                    //         "Equation_Object::build_PADofMap: Number adofs = 0. T-matrix can not be created. MSI probably not build yet. ");

                    MORIS_ASSERT( mFreePdofListFlag,
                            "Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. " );

                    // Get MAX number of pdofs for this equation object
                    uint tNumMyPdofs = mFreePdofList( Ik )( Ij ).size();

                    // MORIS_ASSERT( tNumMyPdofs != 0,
                    //        "Equation_Object::build_PADofMap: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. ");

                    adPADofMap( Ik )( Ij ).resize( tNumMyPdofs );
                    adPADofMapInds( Ik )( Ij ).resize( tNumMyPdofs );

                    // Loop over all pdofs of this equation object
                    for ( uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
                    {
                        auto tPdof = mFreePdofList( Ik )( Ij )( Ii );

                        adPADofMap( Ik )( Ij )( Ii ).resize( tNumUniqueAdofs );
                        adPADofMapInds( Ik )( Ij )( Ii ).resize( tNumUniqueAdofs );

                        // declare a map for pdv-volume fraction derivative.
                        std::unordered_map< moris_id, real > tVolumeFractionDerivMap;

                        // get the numerator and denominator cluster indices for basis extension weights
                        Vector< Vector< moris_index > >  tNumeratorClusterIndices      = ( mFreePdofList( Ik )( Ij )( Ii )->mNumeratorClusterIndices );
                        Vector< Vector< moris_index > >  tDenominatorClusterIndices    = ( mFreePdofList( Ik )( Ij )( Ii )->mDenominatorClusterIndices );
                        Vector< Vector< uint > >         tNumeratorClusterOffsets      = ( mFreePdofList( Ik )( Ij )( Ii )->mNumeratorClusterOffsets );
                        Vector< Vector< uint > >         tDenominatorClusterOffsets    = ( mFreePdofList( Ik )( Ij )( Ii )->mDenominatorClusterOffsets );
                        Vector< Vector< real > >         tBasisWeightsPerLeaderBasis   = ( mFreePdofList( Ik )( Ij )( Ii )->mBasisWeightPerLeaderBasis );

                        // if the size of these above vectors is zero (no extension), fill entries with zero.
                        // fill offsets with zeros as well
                        if ( tNumeratorClusterIndices.size() == 0 )
                        {
                            for ( uint Ib = 0; Ib < tPdof->mAdofIds.numel(); Ib++ )
                            {
                                uint tColumnPos                            = mUniqueAdofMapList( Ik )( Ij )[ tPdof->mAdofIds( Ib, 0 ) ];
                                adPADofMap( Ik )( Ij )(Ii)( tColumnPos ).push_back( 0.0 );
                                adPADofMapInds( Ik )( Ij )(Ii)( tColumnPos ).push_back( gNoIndex );  
                            }
                        }
                        else
                        {
                            // Loop over all adof Ids of this pdof
                            for ( uint Ib = 0; Ib < tPdof->mAdofIds.numel(); Ib++ )
                            {
                                // Getting tPADofMap column entry for the corresponding value
                                uint tColumnPos = mUniqueAdofMapList( Ik )( Ij )[ tPdof->mAdofIds( Ib, 0 ) ];

                                // there are three possibilities here
                                // 1. The SPG was marked for extension, but the basis did not need extension - empty inner vecs in that case
                                // 2. The SPG was marked for extension, and the basis did need extension but no volume fraction averaging needed -> -1 for num and den and offset is size 1
                                // 3. The SPG was marked for extension, and the basis did need extension and volume fraction averaging needed

                                Vector< moris_index > tNumeratorClusterIndicesPerAdof = tNumeratorClusterIndices( Ib );
                                Vector< moris_index > tDenominatorClusterIndicesPerAdof = tDenominatorClusterIndices( Ib );
                                Vector< uint > tNumeratorClusterOffsetsPerAdof = tNumeratorClusterOffsets( Ib );
                                Vector< uint > tDenominatorClusterOffsetsPerAdof = tDenominatorClusterOffsets( Ib );
                                Vector< real > tBasisWeightsPerLeaderBasisPerAdof = tBasisWeightsPerLeaderBasis( Ib );

                                // reset the map to clear the derivatives from the previous numerator denominator pairing
                                tVolumeFractionDerivMap.clear();

                                // to eliminate case 1:
                                if( tNumeratorClusterOffsetsPerAdof.size() >= 1 && tDenominatorClusterOffsetsPerAdof.size() >= 1 )
                                {
                                    // group numerator and denominator clusters into groups signified by the offset
                                    uint tNumeratorClusterOffset = 0;
                                    uint tDenominatorClusterOffset = 0;

                                    // loop over the offset vector to loop over each numerator-denominator pair in the cluster
                                    for ( uint iPair = 0; iPair < tNumeratorClusterOffsetsPerAdof.size(); iPair++ )
                                    {
                                        const auto& tNumeratorClusters = tNumeratorClusterIndicesPerAdof;
                                        const auto& tDenominatorClusters = tDenominatorClusterIndicesPerAdof;
                                        uint tNumCount = tNumeratorClusterOffsetsPerAdof( iPair );
                                        uint tDenCount = tDenominatorClusterOffsetsPerAdof( iPair );

                                        Vector< moris_index > tNumeratorClusterIndicesPerPair;
                                        Vector< moris_index > tDenominatorClusterIndicesPerPair;

                                        tNumeratorClusterIndicesPerPair.insert(
                                                0,
                                                tNumeratorClusters.begin() + tNumeratorClusterOffset,
                                                tNumeratorClusters.begin() + tNumeratorClusterOffset + tNumCount );
                                        tNumeratorClusterOffset += tNumCount;

                                        tDenominatorClusterIndicesPerPair.insert(
                                                0,
                                                tDenominatorClusters.begin() + tDenominatorClusterOffset,
                                                tDenominatorClusters.begin() + tDenominatorClusterOffset + tDenCount );
                                        tDenominatorClusterOffset += tDenCount;

                                        // compute areas and derivative of areas of all clusters in the numerator-denominator pair
                                        Vector< real > tNumClusterAreas( tNumeratorClusterIndicesPerPair.size() );
                                        std::unordered_map< moris_id, real > tNumeratorClusterAreaDerivToPDV;

                                        for ( uint iCluster = 0; iCluster < tNumeratorClusterIndicesPerPair.size(); iCluster++ )
                                        {
                                            // skip if there is only a single spg involved in the volume calculations
                                            if ( tNumeratorClusterIndicesPerPair( iCluster ) == -1 )
                                            {
                                                continue;
                                            }

                                            // compute areas
                                            tNumClusterAreas( iCluster ) = mEquationSet->compute_cluster_measure( tNumeratorClusterIndicesPerPair( iCluster ) );

                                            // compute derivatives with respect to each IG node in the cluster, along each spatial dimension
                                            Matrix< DDRMat > tNumeratorClusterMeasureDerivative =
                                                    mEquationSet->compute_cluster_measure_derivative( tNumeratorClusterIndicesPerPair( iCluster ) );
                                            
                                            // find the list of IG nodes in the cluster
                                            Matrix< DDSMat > tPdvIds = mEquationSet->get_pdvs_in_cluster( tNumeratorClusterIndicesPerPair( iCluster ) );
                                            Matrix< DDSMat > tIsActiveDv = mEquationSet->get_pdvs_activity_in_cluster( tNumeratorClusterIndicesPerPair( iCluster ) );

                                            for (uint iIgNodes = 0 ; iIgNodes < tNumeratorClusterMeasureDerivative.n_rows(); iIgNodes++)

                                            // add the sensitivity of the cluster to the PDV (each PDV represents a spatial coordinate of an IG node)
                                            for (uint iSpatialDim = 0; iSpatialDim < tNumSpatialDims; iSpatialDim++)
                                            {
                                                if ( tPdvIds ( iIgNodes, iSpatialDim ) != -1 && tIsActiveDv( iIgNodes, iSpatialDim ) == 1 )
                                                {
                                                    if ( tNumeratorClusterAreaDerivToPDV.find( static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ) == tNumeratorClusterAreaDerivToPDV.end() )
                                                    {
                                                        tNumeratorClusterAreaDerivToPDV[ static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ] = static_cast< real >( tNumeratorClusterMeasureDerivative( iIgNodes, iSpatialDim ) );
                                                    }
                                                    else
                                                    {
                                                        tNumeratorClusterAreaDerivToPDV[ static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ] += static_cast< real >( tNumeratorClusterMeasureDerivative( iIgNodes, iSpatialDim ) );
                                                    }
                                                        
                                                }
                                            }
                                        }

                                        // compute areas and derivative of areas of all clusters in the numerator-denominator pair
                                        Vector< real > tDenClusterAreas( tDenominatorClusterIndicesPerPair.size() );
                                        std::unordered_map< moris_id, real > tDenominatorClusterAreaDerivToPDV;

                                        for ( uint iCluster = 0; iCluster < tDenominatorClusterIndicesPerPair.size(); iCluster++ )
                                        {
                                            // skip if there is only a single spg involved in the volume calculations
                                            if ( tDenominatorClusterIndicesPerPair( iCluster ) == -1 )
                                            {
                                                continue;
                                            }

                                            // compute areas
                                            tDenClusterAreas( iCluster ) = mEquationSet->compute_cluster_measure( tDenominatorClusterIndicesPerPair( iCluster ) );

                                            // compute derivatives with respect to each IG node in the cluster, along each spatial dimension
                                            Matrix< DDRMat > tDenominatorClusterMeasureDerivative =
                                                    mEquationSet->compute_cluster_measure_derivative( tDenominatorClusterIndicesPerPair( iCluster ) );
                                            
                                            // find the list of IG nodes in the cluster
                                            Matrix< DDSMat > tPdvIds = mEquationSet->get_pdvs_in_cluster( tDenominatorClusterIndicesPerPair( iCluster ) );
                                            Matrix< DDSMat > tIsActiveDv = mEquationSet->get_pdvs_activity_in_cluster( tDenominatorClusterIndicesPerPair( iCluster ) );
                                            

                                            for (uint iIgNodes = 0 ; iIgNodes < tDenominatorClusterMeasureDerivative.n_rows(); iIgNodes++)

                                            // add the sensitivity of the cluster to the PDV (each PDV represents a spatial coordinate of an IG node)
                                            for (uint iSpatialDim = 0; iSpatialDim < tNumSpatialDims; iSpatialDim++)
                                            {
                                                if ( tPdvIds ( iIgNodes, iSpatialDim ) != -1 && tIsActiveDv( iIgNodes, iSpatialDim ) == 1 )
                                                {
                                                    if ( tDenominatorClusterAreaDerivToPDV.find( static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ) == tDenominatorClusterAreaDerivToPDV.end() )
                                                    {
                                                        tDenominatorClusterAreaDerivToPDV[ static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ] = static_cast< real >( tDenominatorClusterMeasureDerivative( iIgNodes, iSpatialDim ) );
                                                    }
                                                    else
                                                    {
                                                        tDenominatorClusterAreaDerivToPDV[ static_cast< moris_id >( tPdvIds( iIgNodes, iSpatialDim ) ) ] += static_cast< real >( tDenominatorClusterMeasureDerivative( iIgNodes, iSpatialDim ) );
                                                    }
                                                        
                                                }
                                            }
                                        }

                                        // now we actually compute the derivative of the ratio with respect to each PDV
                                        for ( auto& tDenominatorClusterAreaDeriv : tDenominatorClusterAreaDerivToPDV )
                                        {
                                            // get the PDV ID
                                            moris_id tPdvId = tDenominatorClusterAreaDeriv.first;

                                            // get the derivative of denominator clusters w.r.t. to this PDV
                                            real tDenDeriv = tDenominatorClusterAreaDeriv.second;

                                            // check if this vertex is also in the numerator IG cells
                                            if ( tNumeratorClusterAreaDerivToPDV.find( tPdvId ) != tNumeratorClusterAreaDerivToPDV.end() )
                                            {
                                                // get the derivative of numerator IG cells w.r.t. to this PDV
                                                real tNumDeriv = tNumeratorClusterAreaDerivToPDV[ tPdvId ];

                                                const real tNumArea =
                                                        std::accumulate( tNumClusterAreas.begin(), tNumClusterAreas.end(), real( 0 ) );
                                                const real tDenArea =
                                                        std::accumulate( tDenClusterAreas.begin(), tDenClusterAreas.end(), real( 0 ) );

                                                // compute the derivative of the volume fraction w.r.t this vertex using the quotient rule
                                                real tVolumeFractionDeriv = ( tNumDeriv * tDenArea - tDenDeriv * tNumArea ) / ( tDenArea * tDenArea );

                                                // compute derivative of b-spline weight this numerator-denominator pairing.
                                                if ( tVolumeFractionDerivMap.find( tPdvId ) == tVolumeFractionDerivMap.end() )
                                                {
                                                    tVolumeFractionDerivMap[ tPdvId ] =  tBasisWeightsPerLeaderBasisPerAdof( iPair ) * tVolumeFractionDeriv;
                                                }
                                                else
                                                {
                                                    tVolumeFractionDerivMap[ tPdvId ] += tBasisWeightsPerLeaderBasisPerAdof( iPair ) * tVolumeFractionDeriv;
                                                }
                                            }
                                            else
                                            {
                                                // No contribution from numerator, so derivative is just the negative of the denominator contribution
                                                const real tNumArea =
                                                        std::accumulate( tNumClusterAreas.begin(), tNumClusterAreas.end(), real( 0 ) );
                                                const real tDenArea =
                                                        std::accumulate( tDenClusterAreas.begin(), tDenClusterAreas.end(), real( 0 ) );

                                                real tVolumeFractionDeriv = ( -tDenDeriv * tNumArea ) / ( tDenArea * tDenArea );

                                                // compute derivative of b-spline weight this numerator-denominator pairing.
                                                if ( tVolumeFractionDerivMap.find( tPdvId ) == tVolumeFractionDerivMap.end() )
                                                {
                                                    tVolumeFractionDerivMap[ tPdvId ] = tBasisWeightsPerLeaderBasisPerAdof( iPair ) * tVolumeFractionDeriv;
                                                }
                                                else
                                                {
                                                    tVolumeFractionDerivMap[ tPdvId ] += tBasisWeightsPerLeaderBasisPerAdof( iPair ) * tVolumeFractionDeriv;
                                                }
                                            }
                                        }
                                    }
                                }

                                // now populate the computed derivatives into the map. If there are zero derivatives (either case 2 or 3) populate the entry with zero.
                                //  populate the output vectors with the IG IDs and their corresponding sensitivities for the current B-Spline basis
                                // if the map is empty (which means either case 2 or case 3 with no other extensions), just fill up with zeros

                                if ( tVolumeFractionDerivMap.empty() )
                                {
                                    adPADofMap( Ik )( Ij )(Ii)( tColumnPos ).push_back( 0.0 );
                                    adPADofMapInds( Ik )( Ij )(Ii)( tColumnPos ).push_back( gNoIndex );
                                }
                                else
                                {
                                    // populate the new map
                                    
                                    for ( const auto& tPair : tVolumeFractionDerivMap )
                                    {
                                        moris_id tPdvId               = tPair.first;
                                        real     tVolumeFractionDeriv = tPair.second;
                                        adPADofMap( Ik )( Ij )(Ii)( tColumnPos ).push_back( tVolumeFractionDeriv );
                                        adPADofMapInds( Ik )( Ij )(Ii)( tColumnPos ).push_back( tPdvId );
                                    }
                                }
                            }
                        }
                        // Loop over the B Spline DoFs and fill zeros in every empty entry (this is because the vectors are
                        // resized as per the total number of unique ADOFs but not every PDoF is associated with every ADOF)
                        for ( uint Ib = 0; Ib < tNumUniqueAdofs; Ib++ )
                        {
                            if ( adPADofMap( Ik )( Ij )(Ii)( Ib ).empty() )
                            {
                                
                                adPADofMap( Ik )( Ij )(Ii)( Ib ).push_back( 0.0 );
                                adPADofMapInds( Ik )( Ij )(Ii)( Ib ).push_back( gNoIndex );
                                
                            }
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_dPADofMap( Vector< Matrix< DDRMat > >& adPADofMap,
                                          Vector< moris_id >& adPADofMapInds )
        {
            /***
             * The goal with this function is to ensure that the large list that was built up in
             * build_PADofMap_list is combined into a single Vector of matrices.
             * The idea is to store it as a matrix corresponding to each PDV, so 
             * the vector loops over the PDVs and the matrices are the dM/dPDV
             */
            Vector< Vector< Vector< Vector< Vector< real > > > > > tdPADofMapList;
            Vector< Vector< Vector< Vector< Vector< moris_id > > > > >    tdPADofMapInds;
            

            // get list of all T-Matrices
            this->build_dPADofMap_list( tdPADofMapList, tdPADofMapInds );

            Vector< enum MSI::Dof_Type > tRequestedDofTypes;

            // get list of requested dof types
            if ( !mEquationSet->mIsStaggered )
            {
                tRequestedDofTypes = mEquationSet->get_requested_dof_types();
            }
            else
            {
                tRequestedDofTypes = mEquationSet->get_secondary_dof_types();
            }

            // initialize column and row counter
            uint tNumColCounter = 0;
            uint tNumRowCounter = 0;
            uint tPdvIdsCounter = 0;

            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                // Loop over all requested dof types and get total number of cols and rows
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get index corresponding to this dof type
                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    // compute the number of PDofs corresponding to a given system and DoF type
                    // important: we assume here that while building the list, the number
                    // of columns is the same for as the vectors have been reshaped to size tNumUniqueAdofs
                    tNumRowCounter += tdPADofMapList( Ii )( tDofTypeIndex ).size();
                    tNumColCounter += mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();


                    // loop over all adofs to find out the total number of PDV IDs associated.
                    for ( uint Ij = 0; Ij < tdPADofMapList( Ii )( tDofTypeIndex ).size(); Ij++ )
                    {
                        for ( uint Il = 0; Il < tdPADofMapList( Ii )( tDofTypeIndex )( Ij ).size(); Il++ )
                        {
                            if ( tdPADofMapInds( Ii )( tDofTypeIndex )( Ij )( Il ).size() > 0 )
                            {
                                tPdvIdsCounter += tdPADofMapInds( Ii )( tDofTypeIndex )( Ij )( Il ).size();
                            }
                        }   
                    }
                }
            }
            // clear previous call data
            adPADofMap.clear();
            adPADofMapInds.clear();

            // initialize
            adPADofMap.resize( tPdvIdsCounter );
            adPADofMapInds.resize( tPdvIdsCounter );

            // re-initialize column and row counter
            uint tNumColCounterNew = 0;
            uint tNumRowCounterNew  = 0;

            std::unordered_map< moris_id, Matrix< DDRMat > > tPdvIdToTmatrixMap;

            // loop to create an unordered map between PDV ID and corresponding T-Matrix derivative
            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get index corresponding to this dof type
                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    // get number of rows and columns
                    uint tNumRows = tdPADofMapList( Ii )( tDofTypeIndex ).size();
                    uint tNumCols = mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();
                    
                    // allocate matrix
                    Matrix< DDRMat > tTmatrix = Matrix< DDRMat >( tNumRowCounter, tNumColCounter, 0.0 );

                    if ( tNumCols * tNumRows > 0 )
                    {
                        // loop over pdofs
                        for ( uint Ij = 0; Ij < tNumRows; Ij++ )
                        {
                            // loop over adofs
                            for ( uint Il = 0; Il < tNumCols; Il++ )
                            {
                                // loop over PDV IDs
                                for ( uint Im = 0; Im < tdPADofMapList( Ii )( tDofTypeIndex )( Ij )( Il ).size(); Im++ )
                                {
                                    // get pdv id
                                    moris_id tPdvId = tdPADofMapInds( Ii )( tDofTypeIndex )( Ij )( Il )( Im );

                                    if ( tPdvId == gNoIndex )
                                    {
                                        continue;
                                    }

                                    // add it to the map if there is nothing there. skip this step if there's already a matrix allocated
                                    if ( tPdvIdToTmatrixMap.find( tPdvId ) == tPdvIdToTmatrixMap.end() )
                                    {
                                        tPdvIdToTmatrixMap[ tPdvId ] = tTmatrix;
                                    }
                                    
                                    // populate the matrix
                                    tPdvIdToTmatrixMap[ tPdvId ]( tNumRowCounterNew + Ij, tNumColCounterNew + Il ) = tdPADofMapList( Ii )( tDofTypeIndex )( Ij )( Il )( Im );
                                        
                                }
                            
                            }
                            
                        }
                    }
                    tNumColCounterNew += tNumCols;
                    tNumRowCounterNew += tNumRows;
                } // end of loop over requested DoF types
            } // end of loop over the requested equation systems: only 1 if bulk element.

            // populate adPADoF
            uint tCounter = 0;
            for ( auto& tEntry : tPdvIdToTmatrixMap )
            {
                adPADofMap( tCounter ) = tEntry.second;
                adPADofMapInds( tCounter ) = tEntry.first;
                tCounter++;
            }
            // reshape the adPADofMap and adPADofMapInds to the correct dimensions
            adPADofMap.resize( tCounter );
            adPADofMapInds.resize( tCounter );
            
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::build_dPADofMap_all( Vector< Matrix< DDRMat > >& adPADofMap,
                                              Vector< moris_id >& adPADofMapInds )
        {
            Vector< Vector< Vector< Vector< Vector< real > > > > > tdPADofMapList;
            Vector< Vector< Vector< Vector< Vector< moris_id > > > > > tdPADofMapInds;

            // Reuse the volume-fraction derivatives and PDV IDs for all systems and dof types.
            this->build_dPADofMap_list( tdPADofMapList, tdPADofMapInds );

            uint tNumMyPdofs = mFreePdofs.size();
            uint tNumUniqueAdofs = mUniqueAdofList.numel();

            MORIS_ASSERT( tNumMyPdofs != 0,
                    "Equation_Object::build_dPADofMap_all: Number pdof types = 0. T-matrix can not be created. MSI probably not build yet. " );

            adPADofMap.clear();
            adPADofMapInds.clear();

            std::unordered_map< moris_id, Matrix< DDRMat > > tPdvIdToTmatrixMap;
            uint tNumRowCounter = 0;

            // The system, dof type and pdof traversal matches create_my_pdof_list.
            for ( uint Ii = 0; Ii < mFreePdofList.size(); Ii++ )
            {
                for ( uint Ij = 0; Ij < mFreePdofList( Ii ).size(); Ij++ )
                {
                    for ( uint Ik = 0; Ik < mFreePdofList( Ii )( Ij ).size(); Ik++ )
                    {
                        auto tPdof = mFreePdofList( Ii )( Ij )( Ik );

                        MORIS_ASSERT( tNumRowCounter < tNumMyPdofs,
                                "Equation_Object::build_dPADofMap_all: Inconsistent number of pdofs." );
                        MORIS_ASSERT( tPdof == mFreePdofs( tNumRowCounter ),
                                "Equation_Object::build_dPADofMap_all: Inconsistent pdof ordering." );

                        for ( uint Ib = 0; Ib < tPdof->mAdofIds.numel(); Ib++ )
                        {
                            moris_id tAdofId = tPdof->mAdofIds( Ib, 0 );
                            uint tLocalColumnPos = mUniqueAdofMapList( Ii )( Ij )[ tAdofId ];
                            uint tColumnPos = mUniqueAdofMap[ tAdofId ];

                            const auto& tDerivatives = tdPADofMapList( Ii )( Ij )( Ik )( tLocalColumnPos );
                            const auto& tPdvIds = tdPADofMapInds( Ii )( Ij )( Ik )( tLocalColumnPos );

                            MORIS_ASSERT( tDerivatives.size() == tPdvIds.size(),
                                    "Equation_Object::build_dPADofMap_all: Derivative values and PDV IDs do not match." );

                            for ( uint Il = 0; Il < tPdvIds.size(); Il++ )
                            {
                                moris_id tPdvId = tPdvIds( Il );

                                if ( tPdvId == gNoIndex )
                                {
                                    continue;
                                }

                                if ( tPdvIdToTmatrixMap.find( tPdvId ) == tPdvIdToTmatrixMap.end() )
                                {
                                    tPdvIdToTmatrixMap[ tPdvId ].set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );
                                }

                                tPdvIdToTmatrixMap[ tPdvId ]( tNumRowCounter, tColumnPos ) = tDerivatives( Il );
                            }
                        }
                        tNumRowCounter++;
                    }
                }
            }

            MORIS_ASSERT( tNumRowCounter == tNumMyPdofs,
                    "Equation_Object::build_dPADofMap_all: Inconsistent number of pdofs." );

            adPADofMap.resize( tPdvIdToTmatrixMap.size() );
            adPADofMapInds.resize( tPdvIdToTmatrixMap.size() );

            uint tCounter = 0;
            for ( const auto& tEntry : tPdvIdToTmatrixMap )
            {
                adPADofMap( tCounter ) = tEntry.second;
                adPADofMapInds( tCounter ) = tEntry.first;
                tCounter++;
            }
        }

        //-------------------------------------------------------------------------------------------------

        moris_index
        Equation_Object::get_node_index( const moris_index aElementLocalNodeIndex ) const
        {
            return mNodeObj( 0 )( aElementLocalNodeIndex )->get_index();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_equation_obj_dof_ids( Matrix< DDSMat >& aEqnObjAdofId )
        {
            Vector< enum MSI::Dof_Type > tRequestedDofTypes;

            // get list of requested dof tpes
            if ( !mEquationSet->mIsStaggered )
            {
                tRequestedDofTypes = mEquationSet->get_requested_dof_types();
            }
            else
            {
                tRequestedDofTypes = mEquationSet->get_secondary_dof_types();
            }

            Dof_Manager* tDofManager = mEquationSet->get_model_solver_interface()->get_dof_manager();

            uint tCounter = 0;
            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    sint tDofTypeIndex =
                            tDofManager->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    tCounter += mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();
                }
            }

            aEqnObjAdofId.set_size( tCounter, 1, -1 );

            tCounter = 0;
            for ( uint Ii = 0; Ii < mNumPdofSystems; Ii++ )
            {
                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    sint tDofTypeIndex = tDofManager->get_pdof_index_for_type( tRequestedDofTypes( Ik ) );

                    uint tNumEntries = mUniqueAdofTypeList( Ii )( tDofTypeIndex ).numel();

                    if ( tNumEntries != 0 )
                    {
                        aEqnObjAdofId( { tCounter, tCounter + tNumEntries - 1 } ) =
                                mUniqueAdofTypeList( Ii )( tDofTypeIndex ).matrix_data();

                        tCounter += tNumEntries;
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_egn_obj_jacobian( Matrix< DDRMat >& aEqnObjMatrix )
        {
            // compute jacobian
            this->compute_jacobian();

            // build T-matrix
            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap_1( tTMatrix );

            // project pdof residual to adof residual
            aEqnObjMatrix = trans( tTMatrix ) * mEquationSet->get_jacobian() * tTMatrix;

            // transpose for sensitivity analysis FIXME move to solver
            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                aEqnObjMatrix = trans( aEqnObjMatrix );
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_equation_obj_residual( Vector< Matrix< DDRMat > >& aEqnObjRHS )
        {
            // compute R for forward analysis or dQIdp for sensitivity analysis
            this->compute_residual();

            // get R or dQIdp values from set
            Vector< Matrix< DDRMat > >& tElementalResidual = mEquationSet->get_residual();

            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                for ( uint Ik = 0; Ik < tElementalResidual.size(); Ik++ )
                {
                    tElementalResidual( Ik ) = -1.0 * tElementalResidual( Ik );
                }

                // this->add_staggered_contribution_to_residual( tElementalResidual );
            }

            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap_1( tTMatrix );

            uint tNumRHS = mEquationSet->mEquationModel->get_num_rhs();

            aEqnObjRHS.resize( tNumRHS );

            // build transpose of Tmatrix
            Matrix< DDRMat > tTMatrixTrans = trans( tTMatrix );

            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {
                aEqnObjRHS( Ik ) = tTMatrixTrans * tElementalResidual( Ik );
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_staggered_equation_obj_residual( Vector< Matrix< DDRMat > >& aEqnObjRHS )
        {
            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                // compute jacobian
                this->compute_jacobian();
            }

            // get R or dQIdp values from set
            Vector< Matrix< DDRMat > >& tElementalResidual = mEquationSet->get_residual();

            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                this->add_staggered_contribution_to_residual( tElementalResidual );
            }

            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap_1( tTMatrix );

            uint tNumRHS = mEquationSet->mEquationModel->get_num_rhs();

            aEqnObjRHS.resize( tNumRHS );

            // build transpose of Tmatrix
            Matrix< DDRMat > tTMatrixTrans = trans( tTMatrix );

            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {

                //MORIS_ASSERT( ( tElementalResidual( Ik ).numel() != 0 ) == ( tTMatrixTrans.numel() != 0 ),
                //        "Equation_Object::get_staggered_equation_obj_residual(), elemental residual vector # %-5i has 0 entries on set %s",
                //        Ik,
                //        mEquationSet->get_set_name().c_str() );

                aEqnObjRHS( Ik ) = tTMatrixTrans * tElementalResidual( Ik );
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_equation_obj_off_diagonal_residual(
                Vector< Matrix< DDRMat > >& aEqnObjRHS )
        {
            // get number of rhs
            uint tNumRHS = mEquationSet->mEquationModel->get_num_rhs();

            // init elemental residual list
            Vector< Matrix< DDRMat > > tElementalResidual( tNumRHS );

            // FIXME this is a hack and will be changed in the next days
            // if sensitivity analysis
            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                // compute jacobian
                this->compute_jacobian();

                // check for zero-size Jacobian
                if ( mEquationSet->get_jacobian().numel() == 0 )
                {
                    aEqnObjRHS.resize( 0 );
                    return;
                }

                // compute previous adjoint values
                this->compute_my_previous_adjoint_values();

                const Vector< enum MSI::Dof_Type >& tRequestedDofTypes = mEquationSet->get_requested_dof_types();

                // get leader dof type list from set
                Vector< Vector< MSI::Dof_Type > >& tLeaderDofTypeList =
                        mEquationSet->get_dof_type_list( mtk::Leader_Follower::LEADER );

                Matrix< DDSMat > tIdentifierMat( (uint)MSI::Dof_Type::END_ENUM, 1, -1 );

                uint tCounter = 0;

                for ( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
                {
                    // get the set index for the requested leader dof type
                    sint tDofIndex = mEquationSet->get_dof_index_for_type(
                            tRequestedDofTypes( Ik ),
                            mtk::Leader_Follower::LEADER );

                    // if the index was set (and is different from -1)
                    if ( tDofIndex != -1 )
                    {
                        tIdentifierMat( tDofIndex ) = 1;

                        uint tEndRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                        tCounter = std::max( tCounter, tEndRow + 1 );
                    }
                }

                Vector< Matrix< DDRMat > > tCoeff( tNumRHS );
                for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                {
                    tCoeff( Ik ).set_size( tCounter, 1, MORIS_REAL_MAX );
                }

                for ( uint Ik = 0; Ik < tLeaderDofTypeList.size(); Ik++ )
                {
                    // get the set index for the requested leader dof type
                    sint tDofIndex = mEquationSet->get_dof_index_for_type(
                            tLeaderDofTypeList( Ik )( 0 ),
                            mtk::Leader_Follower::LEADER );

                    if ( tIdentifierMat( tDofIndex ) == 1 )
                    {
                        // get the pdof values for the ith dof type group
                        Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;

                        this->get_my_pdof_values(
                                mEquationSet->mPreviousAdjointPdofValues,
                                tLeaderDofTypeList( Ik ),
                                tCoeff_Original,
                                mtk::Leader_Follower::LEADER );

                        // loop over the rhs
                        for ( uint Ia = 0; Ia < tNumRHS; Ia++ )
                        {

                            uint tCols = tCoeff_Original( Ia ).size();
                            uint tRows = tCoeff_Original( Ia )( 0 ).numel();

                            uint tStartRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 0 );
                            // uint tEndRow   = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                            for ( uint Ib = 0; Ib < tCols; Ib++ )
                            {
                                for ( uint Ii = 0; Ii < tRows; Ii++ )
                                {
                                    tCoeff( Ia )( tStartRow++ ) = tCoeff_Original( Ia )( Ib )( Ii );
                                }
                            }
                        }
                    }
                }

                // build transpose of Tmatrix
                Matrix< DDRMat > tJacobianTrans = trans( mEquationSet->get_jacobian() );

                // loop over the rhs
                for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                {
                    tElementalResidual( Ik ) = tJacobianTrans * tCoeff( Ik );
                }
            }

            // get the T matrix for eq obj
            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap_1( tTMatrix );

            // init
            aEqnObjRHS.resize( tNumRHS );

            // build transpose of Tmatrix
            Matrix< DDRMat > tTMatrixTrans = trans( tTMatrix );

            // loop over rhs
            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {
                // project
                aEqnObjRHS( Ik ) = tTMatrixTrans * tElementalResidual( Ik );
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::get_egn_obj_jacobian_and_residual(
                Matrix< DDRMat >&         aEqnObjMatrix,
                Vector< Matrix< DDRMat > >& aEqnObjRHS )
        {
            // compute Jacobian and residual
            this->compute_jacobian_and_residual();

            // check for zero-size Jacobian
            // note: if size of Jacobian is zero, also residuals are ignored
            if ( mEquationSet->get_jacobian().numel() == 0 )
            {
                aEqnObjMatrix.set_size( 0, 0 );
                aEqnObjRHS.resize( 0 );
                return;
            }

            // build T-matrix
            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap_1( tTMatrix );

            // get the transpose of the T-matrix (compute once use twice)
            Matrix< DDRMat > tTMatrixTrans = trans( tTMatrix );

            // project pdof residual to adof residual
            aEqnObjMatrix = tTMatrixTrans * mEquationSet->get_jacobian() * tTMatrix;

            // transpose for sensitivity analysis FIXME move to solver
            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                aEqnObjMatrix = trans( aEqnObjMatrix );
            }

            Vector< Matrix< DDRMat > >& tElementalResidual = mEquationSet->get_residual();

            if ( !mEquationSet->mEquationModel->is_forward_analysis() )
            {
                for ( uint Ik = 0; Ik < tElementalResidual.size(); Ik++ )
                {
                    // tElementalResidual( Ik ) = trans( mEquationSet->get_jacobian() ) * mAdjointPdofValues( Ik )- tElementalResidual( Ik );
                    tElementalResidual( Ik ) = -1.0 * tElementalResidual( Ik );
                }
            }

            uint tNumRHS = mEquationSet->mEquationModel->get_num_rhs();

            // resize RHS
            aEqnObjRHS.resize( tNumRHS );

            // multiply RHS with T-matrix
            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
            {
                aEqnObjRHS( Ik ) = tTMatrixTrans * tElementalResidual( Ik );
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::add_staggered_contribution_to_residual( Vector< Matrix< DDRMat > >& aElementResidual )
        {
            Vector< enum MSI::Dof_Type > tAllSecDofTypes = mEquationSet->get_secondary_dof_types();

            if ( tAllSecDofTypes.size() != 0 )
            {
                this->compute_my_adjoint_values();
            }

            for ( auto tSecDofTypes : tAllSecDofTypes )
            {
                Vector< enum MSI::Dof_Type > tRequestedDofTypes = mEquationSet->get_requested_dof_types();

                // combined leader follower index
                sint tSecDofIndex = mEquationSet->get_dof_index_for_type( tSecDofTypes, mtk::Leader_Follower::LEADER );

                if ( tSecDofIndex != -1 )
                {
                    for ( auto tDofTypes : tRequestedDofTypes )
                    {
                        sint tDofIndex = mEquationSet->get_dof_index_for_type( tDofTypes, mtk::Leader_Follower::LEADER );

                        if ( tDofIndex != -1 )
                        {
                            uint tStartRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 0 );
                            uint tEndRow   = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                            uint tStartCol = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                            uint tEndCol   = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                            // build transpose of Jacobian
                            Matrix< DDRMat > tJacTrans =
                                    trans( mEquationSet->get_jacobian()(
                                            { tStartRow, tEndRow },
                                            { tStartCol, tEndCol } ) );

                            // separate leader follower index
                            sint tseparateSecDofIndex = mEquationSet->get_dof_index_for_type_1( tDofTypes, mtk::Leader_Follower::LEADER );

                            // get the ith dof type group
                            const Vector< MSI::Dof_Type >& tDofTypeGroup = mEquationSet->mLeaderDofTypes( tseparateSecDofIndex );

                            // get the pdof values for the ith dof type group
                            Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;

                            this->get_my_pdof_values( mEquationSet->mAdjointPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Leader_Follower::LEADER );

                            uint                     tNumRHS = tCoeff_Original.size();
                            Vector< Matrix< DDRMat > > tCoeff( tNumRHS );

                            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                            {
                                // reshape tCoeffs into the order the cluster expects them
                                this->reshape_pdof_values( tCoeff_Original( Ik ), tCoeff( Ik ) );

                                Vector< Matrix< DDRMat > > tCoeff1( tNumRHS );

                                tCoeff1( Ik ).set_size( tCoeff( Ik ).numel(), 1, 0.0 );

                                // fixme get rid of for loop
                                uint tCounter = 0;
                                for ( uint Ia = 0; Ia < tCoeff( Ik ).n_cols(); Ia++ )
                                {
                                    for ( uint Ii = 0; Ii < tCoeff( Ik ).n_rows(); Ii++ )
                                    {
                                        tCoeff1( Ik )( tCounter++ ) = tCoeff( Ik )( Ii, Ia );
                                    }
                                }

                                aElementResidual( Ik )( { tStartCol, tEndCol } ) += tJacTrans * tCoeff1( Ik );
                            }
                        }

                        tDofIndex = mEquationSet->get_dof_index_for_type( tDofTypes, mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex != -1 )
                        {
                            uint tStartRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 0 );
                            uint tEndRow   = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                            uint tStartCol = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                            uint tEndCol   = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                            // build transpose of Jacobian
                            Matrix< DDRMat > tJacTrans =
                                    trans( mEquationSet->get_jacobian()(
                                            { tStartRow, tEndRow },
                                            { tStartCol, tEndCol } ) );

                            // separate leader follower index
                            sint tSeparateSecDofIndex = mEquationSet->get_dof_index_for_type_1( tDofTypes, mtk::Leader_Follower::FOLLOWER );

                            // get the ith dof type group
                            const Vector< MSI::Dof_Type >& tDofTypeGroup = mEquationSet->mFollowerDofTypes( tSeparateSecDofIndex );

                            // get the pdof values for the ith dof type group
                            Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;

                            this->get_my_pdof_values( mEquationSet->mAdjointPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Leader_Follower::FOLLOWER );

                            uint                     tNumRHS = tCoeff_Original.size();
                            Vector< Matrix< DDRMat > > tCoeff( tNumRHS );

                            // print(mEquationSet->get_jacobian(),"staggerd jac");

                            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                            {
                                // reshape tCoeffs into the order the cluster expects them
                                this->reshape_pdof_values( tCoeff_Original( Ik ), tCoeff( Ik ) );

                                Vector< Matrix< DDRMat > > tCoeff1( tNumRHS );

                                tCoeff1( Ik ).set_size( tCoeff( Ik ).numel(), 1, 0.0 );

                                // fixme get rid of for loop
                                uint tCounter = 0;
                                for ( uint Ia = 0; Ia < tCoeff( Ik ).n_cols(); Ia++ )
                                {
                                    for ( uint Ii = 0; Ii < tCoeff( Ik ).n_rows(); Ii++ )
                                    {
                                        tCoeff1( Ik )( tCounter++ ) = tCoeff( Ik )( Ii, Ia );
                                    }
                                }

                                aElementResidual( Ik )( { tStartCol, tEndCol } ) += tJacTrans * tCoeff1( Ik );
                            }
                        }
                    }
                }

                // combined leader follower index
                tSecDofIndex = mEquationSet->get_dof_index_for_type( tSecDofTypes, mtk::Leader_Follower::FOLLOWER );

                if ( tSecDofIndex != -1 )
                {
                    for ( auto tDofTypes : tRequestedDofTypes )
                    {
                        sint tDofIndex = mEquationSet->get_dof_index_for_type( tDofTypes, mtk::Leader_Follower::LEADER );

                        if ( tDofIndex != -1 )
                        {
                            uint tStartRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 0 );
                            uint tEndRow   = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                            uint tStartCol = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                            uint tEndCol   = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                            // build transpose of Jacobian
                            Matrix< DDRMat > tJacTrans =
                                    trans( mEquationSet->get_jacobian()(
                                            { tStartRow, tEndRow },
                                            { tStartCol, tEndCol } ) );

                            // separate leader follower index
                            sint tseparateSecDofIndex = mEquationSet->get_dof_index_for_type_1( tDofTypes, mtk::Leader_Follower::LEADER );

                            // get the ith dof type group
                            const Vector< MSI::Dof_Type >& tDofTypeGroup = mEquationSet->mLeaderDofTypes( tseparateSecDofIndex );

                            // get the pdof values for the ith dof type group
                            Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;

                            this->get_my_pdof_values( mEquationSet->mAdjointPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Leader_Follower::LEADER );

                            uint                     tNumRHS = tCoeff_Original.size();
                            Vector< Matrix< DDRMat > > tCoeff( tNumRHS );

                            // print(mEquationSet->get_jacobian(),"staggerd jac");

                            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                            {
                                // reshape tCoeffs into the order the cluster expects them
                                this->reshape_pdof_values( tCoeff_Original( Ik ), tCoeff( Ik ) );

                                Vector< Matrix< DDRMat > > tCoeff1( tNumRHS );

                                tCoeff1( Ik ).set_size( tCoeff( Ik ).numel(), 1, 0.0 );

                                // fixme get rid of for loop
                                uint tCounter = 0;
                                for ( uint Ia = 0; Ia < tCoeff( Ik ).n_cols(); Ia++ )
                                {
                                    for ( uint Ii = 0; Ii < tCoeff( Ik ).n_rows(); Ii++ )
                                    {
                                        tCoeff1( Ik )( tCounter++ ) = tCoeff( Ik )( Ii, Ia );
                                    }
                                }

                                aElementResidual( Ik )( { tStartCol, tEndCol } ) += tJacTrans * tCoeff1( Ik );
                            }
                        }

                        tDofIndex = mEquationSet->get_dof_index_for_type( tDofTypes, mtk::Leader_Follower::FOLLOWER );

                        if ( tDofIndex != -1 )
                        {
                            uint tStartRow = mEquationSet->mResDofAssemblyMap( tDofIndex )( 0 );
                            uint tEndRow   = mEquationSet->mResDofAssemblyMap( tDofIndex )( 1 );

                            uint tStartCol = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 0 );
                            uint tEndCol   = mEquationSet->mJacDofAssemblyMap( tDofIndex )( tSecDofIndex, 1 );

                            // build transpose of Jacobian
                            Matrix< DDRMat > tJacTrans =
                                    trans( mEquationSet->get_jacobian()(
                                            { tStartRow, tEndRow },
                                            { tStartCol, tEndCol } ) );

                            // separate leader follower index
                            sint tseparateSecDofIndex = mEquationSet->get_dof_index_for_type_1( tDofTypes, mtk::Leader_Follower::FOLLOWER );

                            // get the ith dof type group
                            const Vector< MSI::Dof_Type >& tDofTypeGroup = mEquationSet->mFollowerDofTypes( tseparateSecDofIndex );

                            // get the pdof values for the ith dof type group
                            Vector< Vector< Matrix< DDRMat > > > tCoeff_Original;

                            this->get_my_pdof_values( mEquationSet->mAdjointPdofValues, tDofTypeGroup, tCoeff_Original, mtk::Leader_Follower::FOLLOWER );

                            uint                     tNumRHS = tCoeff_Original.size();
                            Vector< Matrix< DDRMat > > tCoeff( tNumRHS );

                            // print(mEquationSet->get_jacobian(),"staggerd jac");

                            for ( uint Ik = 0; Ik < tNumRHS; Ik++ )
                            {
                                // reshape tCoeffs into the order the cluster expects them
                                this->reshape_pdof_values( tCoeff_Original( Ik ), tCoeff( Ik ) );

                                Vector< Matrix< DDRMat > > tCoeff1( tNumRHS );

                                tCoeff1( Ik ).set_size( tCoeff( Ik ).numel(), 1, 0.0 );

                                // fixme get rid of for loop
                                uint tCounter = 0;
                                for ( uint Ia = 0; Ia < tCoeff( Ik ).n_cols(); Ia++ )
                                {
                                    for ( uint Ii = 0; Ii < tCoeff( Ik ).n_rows(); Ii++ )
                                    {
                                        tCoeff1( Ik )( tCounter++ ) = tCoeff( Ik )( Ii, Ia );
                                    }
                                }

                                aElementResidual( Ik )( { tStartCol, tEndCol } ) += tJacTrans * tCoeff1( Ik );
                            }
                        }
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_pdof_values()
        {
            Matrix< DDRMat > tTMatrix;

            // build T-matrix
            this->build_PADofMap( tTMatrix );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel->get_solution_vector()->    //
                    extract_my_values(
                            tTMatrix.n_cols(),
                            mUniqueAdofList,
                            0,
                            tMyValues );

            if ( mEquationSet->mPdofValues.size() != tMyValues.size() )
            {
                mEquationSet->mPdofValues.resize( tMyValues.size() );
            }

            // multiply t_matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mPdofValues( Ik ) = tTMatrix * tMyValues( Ik );
            }

            this->set_vector_entry_number_of_pdof();    // FIXME should not be in MSI. Should be in FEM
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_previous_pdof_values()
        {
            Matrix< DDRMat > tTMatrix;

            // build T-matrix
            this->build_PADofMap( tTMatrix );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_previous_solution_vector()
                    ->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

            if ( mEquationSet->mPreviousPdofValues.size() != tMyValues.size() )
            {
                mEquationSet->mPreviousPdofValues.resize( tMyValues.size() );
            }

            // multiply t_matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mPreviousPdofValues( Ik ) = tTMatrix * tMyValues( Ik );
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_eigen_vector_values()
        {
            Matrix< DDRMat > tTMatrix;

            // build T-matrix
            this->build_PADofMap( tTMatrix );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_eigen_solution_vector()
                    ->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

            if ( mEquationSet->mEigenVectorPdofValues.size() != tMyValues.size() )
            {
                mEquationSet->mEigenVectorPdofValues.resize( tMyValues.size() );
            }

            // multiply t_matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mEigenVectorPdofValues( Ik ) = tTMatrix * tMyValues( Ik );
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_adjoint_values()
        {
            Matrix< DDRMat > tTMatrix;

            // build T-matrix
            this->build_PADofMap( tTMatrix );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_adjoint_solution_vector()
                    ->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

            if ( mEquationSet->mAdjointPdofValues.size() != tMyValues.size() )
            {
                mEquationSet->mAdjointPdofValues.resize( tMyValues.size() );
            }

            // multiply t_matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mAdjointPdofValues( Ik ) = tTMatrix * tMyValues( Ik );
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_dPADofMap_adjoint_values()
        {
            Vector< Matrix< DDRMat > > tDTMatrix;
            Matrix< DDRMat > tTMatrix;
            Vector< moris_id > tDTMatrixIndices;

            // get T matrix
            this->build_PADofMap( tTMatrix );

            // get derivative of T Matrix
            this->build_dPADofMap_all( tDTMatrix, tDTMatrixIndices );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_adjoint_solution_vector()
                    ->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

            mEquationSet->mdPADofAdjointPdofValues.resize( tMyValues.size() );
            mEquationSet->mdPADofAdjointPdofIndices.resize( tMyValues.size() );

            // multiply derivative of T-matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                // loop over pdvs
                mEquationSet->mdPADofAdjointPdofValues( Ik ).resize( tDTMatrixIndices.size() );
                mEquationSet->mdPADofAdjointPdofIndices( Ik ) = tDTMatrixIndices;

                for ( uint Ij = 0; Ij < tDTMatrixIndices.size(); Ij++ )
                {
                    mEquationSet->mdPADofAdjointPdofValues( Ik )( Ij ) = tDTMatrix( Ij ) * tMyValues( Ik );
                }
                
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_previous_adjoint_values()
        {
            Matrix< DDRMat > tTMatrix;

            // build T-matrix
            this->build_PADofMap( tTMatrix );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel->get_previous_adjoint_solution_vector()->extract_my_values( tTMatrix.n_cols(), mUniqueAdofList, 0, tMyValues );

            if ( mEquationSet->mPreviousAdjointPdofValues.size() != tMyValues.size() )
            {
                mEquationSet->mPreviousAdjointPdofValues.resize( tMyValues.size() );
            }

            // multiply t_matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mPreviousAdjointPdofValues( Ik ) = tTMatrix * tMyValues( Ik );
            }

            this->set_vector_entry_number_of_pdof();    // FIXME should not be in MSI. Should be in FEM
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_dPADofMap_pdof_values()
        {
            Vector< Matrix< DDRMat > > tDTMatrix;
            Vector< moris_id > tDTMatrixIndices;

            // get derivative of T-matrix
            this->build_dPADofMap_all( tDTMatrix, tDTMatrixIndices );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_solution_vector()
                    ->extract_my_values( mUniqueAdofList.numel(), mUniqueAdofList, 0, tMyValues );

            mEquationSet->mdPADofPdofValues.resize( tMyValues.size() );
            mEquationSet->mdPADofPdofIndices.resize( tMyValues.size() );

            // multiply derivative of T-matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mdPADofPdofValues( Ik ).resize( tDTMatrixIndices.size() );
                mEquationSet->mdPADofPdofIndices( Ik ) = tDTMatrixIndices;

                for ( uint Ij = 0; Ij < tDTMatrixIndices.size(); Ij++ )
                {
                    mEquationSet->mdPADofPdofValues( Ik )( Ij ) = tDTMatrix( Ij ) * tMyValues( Ik );
                }
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_previous_dPADofMap_pdof_values()
        {
            Vector< Matrix< DDRMat > > tDTMatrix;
            Vector< moris_id > tDTMatrixIndices;

            // get derivative of T-matrix
            this->build_dPADofMap_all( tDTMatrix, tDTMatrixIndices );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_previous_solution_vector()
                    ->extract_my_values( mUniqueAdofList.numel(), mUniqueAdofList, 0, tMyValues );

            mEquationSet->mdPADofPreviousPdofValues.resize( tMyValues.size() );
            mEquationSet->mdPADofPreviousPdofIndices.resize( tMyValues.size() );

            // multiply derivative of T-matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mdPADofPreviousPdofValues( Ik ).resize( tDTMatrixIndices.size() );
                mEquationSet->mdPADofPreviousPdofIndices( Ik ) = tDTMatrixIndices;

                for ( uint Ij = 0; Ij < tDTMatrixIndices.size(); Ij++ )
                {
                    mEquationSet->mdPADofPreviousPdofValues( Ik )( Ij ) = tDTMatrix( Ij ) * tMyValues( Ik );
                }
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_dPADofMap_eigen_vector_values()
        {
            Vector< Matrix< DDRMat > > tDTMatrix;
            Vector< moris_id > tDTMatrixIndices;

            // get derivative of T-matrix
            this->build_dPADofMap_all( tDTMatrix, tDTMatrixIndices );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_eigen_solution_vector()
                    ->extract_my_values( mUniqueAdofList.numel(), mUniqueAdofList, 0, tMyValues );

            mEquationSet->mdPADofEigenVectorPdofValues.resize( tMyValues.size() );
            mEquationSet->mdPADofEigenVectorPdofIndices.resize( tMyValues.size() );

            // multiply derivative of T-matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mdPADofEigenVectorPdofValues( Ik ).resize( tDTMatrixIndices.size() );
                mEquationSet->mdPADofEigenVectorPdofIndices( Ik ) = tDTMatrixIndices;

                for ( uint Ij = 0; Ij < tDTMatrixIndices.size(); Ij++ )
                {
                    mEquationSet->mdPADofEigenVectorPdofValues( Ik )( Ij ) = tDTMatrix( Ij ) * tMyValues( Ik );
                }
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::compute_my_previous_dPADofMap_adjoint_values()
        {
            Vector< Matrix< DDRMat > > tDTMatrix;
            Vector< moris_id > tDTMatrixIndices;

            // get derivative of T-matrix
            this->build_dPADofMap_all( tDTMatrix, tDTMatrixIndices );

            Vector< Matrix< DDRMat > > tMyValues;

            // Extract this equation objects adof values from solution vector
            mEquationSet->mEquationModel
                    ->get_previous_adjoint_solution_vector()
                    ->extract_my_values( mUniqueAdofList.numel(), mUniqueAdofList, 0, tMyValues );

            mEquationSet->mdPADofPreviousAdjointPdofValues.resize( tMyValues.size() );
            mEquationSet->mdPADofPreviousAdjointPdofIndices.resize( tMyValues.size() );

            // multiply derivative of T-matrix with adof values to get pdof values
            for ( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            {
                mEquationSet->mdPADofPreviousAdjointPdofValues( Ik ).resize( tDTMatrixIndices.size() );
                mEquationSet->mdPADofPreviousAdjointPdofIndices( Ik ) = tDTMatrixIndices;

                for ( uint Ij = 0; Ij < tDTMatrixIndices.size(); Ij++ )
                {
                    mEquationSet->mdPADofPreviousAdjointPdofValues( Ik )( Ij ) = tDTMatrix( Ij ) * tMyValues( Ik );
                }
            }

            // FIXME should not be in MSI. Should be in FEM
            this->set_vector_entry_number_of_pdof();
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::set_time( Matrix< DDRMat >& aTime )
        {
            return mEquationSet->mEquationModel->set_time( aTime );
        }

        //-------------------------------------------------------------------------------------------------

        Matrix< DDRMat >&
        Equation_Object::get_time()
        {
            return mEquationSet->mEquationModel->get_time();
        }

        //-------------------------------------------------------------------------------------------------

        Matrix< DDRMat >&
        Equation_Object::get_previous_time()
        {
            return mEquationSet->mEquationModel->get_previous_time();
        }

        //-------------------------------------------------------------------------------------------------

        //    void Equation_Object::get_my_pdof_values( const Vector< enum Dof_Type > & aRequestedDofTypes,
        //                                                    Cell< Matrix< DDRMat > >     & aRequestedPdofValues )
        //    {
        //        // Initialize list which contains the maximal number of time levels per dof type
        //        Matrix< DDSMat > tTimeLevelsPerDofType( aRequestedDofTypes.size(), 1, -1 );
        //
        //        sint tCounter = 0;
        //
        //        // Loop over requested dof types
        //        for ( uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        //        {
        //            // Loop over all elemental pdof hosts
        //            for ( uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
        //            {
        //                // Get dof type index
        //                sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()
        //                                                                 ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );
        //
        //                MORIS_ERROR( mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) !=0,
        //                        "Equation_Object::get_my_pdof_values: talk with Mathias about this");                         //FIXME delete this error after a closer look
        //
        //                // get number of time levels for this dof type
        //                sint tNumTimeLevels = mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex );
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
        //        sint tCounter_2 = 0;
        //
        //        // Loop over requested dof types
        //        for ( uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
        //        {
        //            // Get maximal Number of time levels on this pdof type
        //            sint tMaxTimeLevelsOnDofType = tTimeLevelsPerDofType( Ii, 0 );
        //
        //            // Loop over this pdofs time levels
        //            for ( sint Ia = 0; Ia < tMaxTimeLevelsOnDofType; Ia++ )
        //            {
        //                // Loop over all elemental pdof hosts
        //                for ( uint Ik = 0; Ik < mMyPdofHosts.size(); Ik++ )
        //                {
        //                    // Get dof type index
        //                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()
        //                                                                     ->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );
        //
        //                    // Check if number if time levels on this dof type is smaller than maximal number of time levels on dof type
        //                    if ( (sint)mMyPdofHosts( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) == tMaxTimeLevelsOnDofType )
        //                    {
        //                        // get pointer list all time pdofs on this pdof type
        //                        Vector< Pdof* > tPdofTimeList = mMyPdofHosts( Ik )->get_pdof_time_list( tDofTypeIndex );
        //
        //                        // get entry number of this pdof in the elemental pdof value vector
        //                        uint mElementalSolVecEntry = tPdofTimeList( Ia )->mElementalSolVecEntry;
        //
        //                        // Put this pdof value into the requested pdof vector
        //                        aRequestedPdofValues( tCounter_2++, 0 ) = mPdofValues( mElementalSolVecEntry , 0 );
        //                    }
        //                }
        //            }
        //        }
        //    }

        void
        Equation_Object::get_my_pdof_values(
                const Vector< Matrix< DDRMat > >&   aPdofValues,
                const Vector< enum Dof_Type >&      aRequestedDofTypes,
                Vector< Vector< Matrix< DDRMat > > >& aRequestedPdofValues,
                const mtk::Leader_Follower                  aIsLeader )
        {
            // check that leader or follower
            MORIS_ERROR( ( aIsLeader == mtk::Leader_Follower::LEADER ) || ( aIsLeader == mtk::Leader_Follower::FOLLOWER ),
                    "Equation_Object::get_my_pdof_values - can only be LEADER or FOLLOWER" );

            // uint for leader/follower
            uint tIsLeader = static_cast< uint >( aIsLeader );

            // Initialize list which contains the maximal number of time levels per dof type
            Matrix< DDSMat > tTimeLevelsPerDofType( aRequestedDofTypes.size(), 1, -1 );

            uint tNumVectors = aPdofValues.size();

            // set size for number of solution vectors
            aRequestedPdofValues.resize( tNumVectors );

            for ( uint Ik = 0; Ik < tNumVectors; Ik++ )
            {
                aRequestedPdofValues( Ik ).resize( aRequestedDofTypes.size() );
            }

            sint tCounter = 0;

            // Loop over requested dof types
            for ( uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
            {
                tCounter = 0;

                // Loop over all elemental pdof hosts
                for ( uint Ik = 0; Ik < mMyPdofHosts( tIsLeader ).size(); Ik++ )
                {
                    // Get dof type index
                    sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                    // FIXME delete this error after a closer look
                    MORIS_ASSERT( mMyPdofHosts( tIsLeader )( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) != 0,
                            "Equation_Object::get_my_pdof_values: talk with Mathias about this" );

                    // get number of time levels for this dof type
                    sint tNumTimeLevels = mMyPdofHosts( tIsLeader )( Ik )->get_num_time_levels_of_type( tDofTypeIndex );
                    tCounter            = tCounter + tNumTimeLevels;

                    // Add maximal value of time levels to list
                    tTimeLevelsPerDofType( Ii, 0 ) = std::max( tTimeLevelsPerDofType( Ii, 0 ), tNumTimeLevels );
                }
                MORIS_ASSERT( tTimeLevelsPerDofType( Ii, 0 ) > -1,
                        "Equation_Object::get_my_pdof_values: no time levels exist on this dof type on element %-5i",
                        mEqnObjInd );

                // set size for all solution vectors
                for ( uint Ia = 0; Ia < tNumVectors; Ia++ )
                {
                    // Set size matrix for requested pdof values
                    aRequestedPdofValues( Ia )( Ii ).resize( tCounter, 1 );
                }
            }

            sint tCounter_2 = 0;

            // Loop over requested dof types
            for ( uint Ii = 0; Ii < aRequestedDofTypes.size(); Ii++ )
            {
                tCounter_2 = 0;
                // Get maximal Number of time levels on this pdof type
                sint tMaxTimeLevelsOnDofType = tTimeLevelsPerDofType( Ii, 0 );

                // Loop over this pdofs time levels
                for ( sint Ia = 0; Ia < tMaxTimeLevelsOnDofType; Ia++ )
                {
                    // Loop over all elemental pdof hosts
                    for ( uint Ik = 0; Ik < mMyPdofHosts( tIsLeader ).size(); Ik++ )
                    {
                        // Get dof type index
                        sint tDofTypeIndex = mEquationSet->get_model_solver_interface()->get_dof_manager()->get_pdof_index_for_type( aRequestedDofTypes( Ii ) );

                        // Check if number of time levels on this dof type is smaller than maximal number of time levels on dof type
                        if ( (sint)mMyPdofHosts( tIsLeader )( Ik )->get_num_time_levels_of_type( tDofTypeIndex ) == tMaxTimeLevelsOnDofType )
                        {
                            // get pointer list all time pdofs on this pdof type
                            Vector< Pdof* > tPdofTimeList = mMyPdofHosts( tIsLeader )( Ik )->get_pdof_time_list( tDofTypeIndex );

                            // get entry number of this pdof in the elemental pdof value vector
                            uint tElementalSolVecEntry = tPdofTimeList( Ia )->mElementalSolVecEntry;

                            for ( uint Ib = 0; Ib < tNumVectors; Ib++ )
                            {
                                // Put this pdof value into the requested pdof vector
                                aRequestedPdofValues( Ib )( Ii )( tCounter_2, 0 ) = aPdofValues( Ib )( tElementalSolVecEntry, 0 );
                            }
                            tCounter_2++;
                        }
                    }
                }
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::reshape_pdof_values(
                const Vector< Matrix< DDRMat > >& aPdofValues,
                Matrix< DDRMat >&               aReshapedPdofValues )
        {
            MORIS_ASSERT( aPdofValues.size() != 0,
                    "Equation_Object::reshape_pdof_values(), pdof value vector is empty" );

            uint tCols = aPdofValues.size();
            uint tRows = aPdofValues( 0 ).numel();

            aReshapedPdofValues.set_size( tRows, tCols );

            for ( uint Ik = 0; Ik < tCols; Ik++ )
            {
                aReshapedPdofValues( { 0, tRows - 1 }, { Ik, Ik } ) =
                        aPdofValues( Ik ).matrix_data();
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::reshape_pdof_values_vector(
                const Vector< Matrix< DDRMat > >& aPdofValues,
                Matrix< DDRMat >&               aReshapedPdofValues )
        {
            MORIS_ASSERT( aPdofValues.size() != 0,
                    "Equation_Object::reshape_pdof_values(), pdof value vector is empty" );

            uint tCols = aPdofValues.size();
            uint tRows = aPdofValues( 0 ).numel();

            aReshapedPdofValues.set_size( tRows * tCols, 1 );

            for ( uint Ik = 0; Ik < tCols; Ik++ )
            {
                aReshapedPdofValues( { Ik * tRows, ( Ik + 1 ) * tRows - 1 } ) =
                        aPdofValues( Ik ).matrix_data();
            }
        }

        //-------------------------------------------------------------------------------------------------

        void
        Equation_Object::set_vector_entry_number_of_pdof()
        {
            uint tNumMyPdofs = mFreePdofs.size();
            // Loop over all pdofs of this element
            for ( uint Ik = 0; Ik < tNumMyPdofs; Ik++ )
            {
                mFreePdofs( Ik )->mElementalSolVecEntry = Ik;
            }
        }
    }    // namespace MSI
}    // namespace moris
