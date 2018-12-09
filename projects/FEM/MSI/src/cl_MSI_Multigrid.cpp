/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Multigrid.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

#include "fn_print.hpp"

namespace moris
{
    namespace MSI
    {
    Multigrid::Multigrid( moris::MSI::Model_Solver_Interface * aModelSolverInterface,
                          moris::mtk::Mesh                   * aMesh ) : mMesh( aMesh ),
                                                                         mModelSolverInterface( aModelSolverInterface )
    {
        mMultigridLevels = 2;   //FIXME Input

        mNumDofsRemain.set_size( mMultigridLevels, 1 );
    }

//------------------------------------------------------------------------------------------------------------------------------------

    void Multigrid::multigrid_initialize()
    {
        moris::Cell < Adof * > tOwnedAdofList = mModelSolverInterface->get_dof_manager()->get_owned_adofs();

        moris::uint tNumOwnedAdofs = tOwnedAdofList.size();

        mListAdofExtIndMap         .resize( mMultigridLevels + 1 );
        mListAdofTypeTimeIdentifier.resize( mMultigridLevels + 1 );

        mListAdofExtIndMap( 0 )         .set_size( tNumOwnedAdofs, 1 );
        mListAdofTypeTimeIdentifier( 0 ).set_size( tNumOwnedAdofs, 1 );

        moris::uint Ik = 0;
        for ( Adof* tAdof : tOwnedAdofList )
        {
            mListAdofExtIndMap( 0 )( Ik, 0) = tAdof->get_adof_external_ind();

            mListAdofTypeTimeIdentifier( 0 )( Ik++, 0) = tAdof->get_adof_type_time_identifier();
        }

        mMaxDofTypes = mListAdofTypeTimeIdentifier( 0 ).max()+1;

        MORIS_ASSERT( mListAdofTypeTimeIdentifier( 0 ).min() != -1, "moris::MSI::Multigrid: Type/time identifier not specified");

        this->create_multigrid_level_dof_ordering();

        this->create_multigrid_maps();
    }

//------------------------------------------------------------------------------------------------------------------------------------

    void Multigrid::create_multigrid_level_dof_ordering()
    {
        // Gets the maximal mesh level
        moris::uint tMaxMeshLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )->get_max_level();

        // Loop over all multigrid levels
        for ( moris::sint Ik = 0; Ik < mMultigridLevels; Ik++ )
        {
            moris::uint tNumDofsOnLevel = mListAdofExtIndMap( Ik ).length();
            moris::uint tMaxIndex = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )
                                                             ->get_number_of_indexed_basis();

            moris::Cell < Matrix < DDUMat> > tCoarseDofExist( mMaxDofTypes );
            for ( moris::sint Ib = 0; Ib < mMaxDofTypes; Ib++ )
            {
                tCoarseDofExist( Ib ).set_size( tMaxIndex, 1, 0 );
            }

            moris::uint tCounter = 0;
            moris::uint tCounterTooFine = 0;

            mListAdofExtIndMap( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );
            mListAdofTypeTimeIdentifier( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );

            Matrix< DDSMat > tEntryOfTooFineDofs( tNumDofsOnLevel, 1, -1);

            // Loop over all dofs on this level
            for ( moris::uint Ii = 0; Ii < tNumDofsOnLevel; Ii++ )
            {
                // Ask mesh for the level of this mesh index
                moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )
                                                                 ->get_basis_by_index( mListAdofExtIndMap( Ik )( Ii, 0 ) )
                                                                 ->get_level();

                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
                if( tDofLevel < tMaxMeshLevel - Ik )
                {
                    mListAdofExtIndMap( Ik + 1 )( tCounter ) = mListAdofExtIndMap( Ik )( Ii, 0 );

                    mListAdofTypeTimeIdentifier( Ik + 1 )( tCounter++ ) = mListAdofTypeTimeIdentifier( Ik )( Ii, 0 );

                    tCoarseDofExist( mListAdofTypeTimeIdentifier( Ik )( Ii, 0 ) )( mListAdofExtIndMap( Ik )( Ii, 0 ), 0 ) = 1;
                }
                else
                {
                    tEntryOfTooFineDofs( tCounterTooFine++, 0 ) = Ii;
                }
            }

            mNumDofsRemain( Ik, 0 ) = tCounter;

            tEntryOfTooFineDofs.resize( tCounterTooFine, 1 );

            // Loop over all refined dofs on this level
            for ( moris::uint Ij = 0; Ij < tCounterTooFine; Ij++ )
            {
                moris::hmr::Basis * tBasis = mMesh->get_HMR_database()
                                                  ->get_bspline_mesh_by_index( gAdofOrderHack )
                                                  ->get_basis_by_index( mListAdofExtIndMap( Ik )( tEntryOfTooFineDofs( Ij, 0 ), 0 ) );

                // Get type/time identifier for this dof
                moris::uint tTypeIdentifier = mListAdofTypeTimeIdentifier( Ik )( tEntryOfTooFineDofs( Ij, 0 ), 0 );

                moris:: uint tNumCoarseDofs = tBasis ->get_number_of_parents() ;

                for ( moris::uint Ia = 0; Ia < tNumCoarseDofs; Ia++ )
                {
                    moris:: uint tCoarseDofIndex = tBasis->get_parent( Ia )->get_index();

                    if ( tCoarseDofExist( tTypeIdentifier )( tCoarseDofIndex , 0 ) == 0 )
                    {
                        mListAdofExtIndMap( Ik + 1 )( tCounter ) = tCoarseDofIndex;

                        mListAdofTypeTimeIdentifier( Ik + 1 )( tCounter++ ) = tTypeIdentifier;

                        tCoarseDofExist( tTypeIdentifier )( tCoarseDofIndex, 0 ) = 1;
                    }
                }
            }
            mListAdofExtIndMap( Ik + 1 ).resize( tCounter, 1 );

            mListAdofTypeTimeIdentifier( Ik + 1 ).resize( tCounter, 1 );
        }
    }

    //------------------------------------------------------------------------------------------------------------------------------------

    void Multigrid::create_multigrid_maps()
    {
        moris::uint tMaxIndex = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )
                                                         ->get_number_of_indexed_basis();

        mMultigridMap.resize( mMultigridLevels + 1 );

        // Loop to set size of list of list of mats
        for ( moris::uint Ik = 0; Ik < mMultigridMap.size(); Ik++ )
        {
            mMultigridMap( Ik ).resize( mMaxDofTypes );

            for ( moris::uint Ii = 0; Ii < mMultigridMap( Ik ).size(); Ii++ )
            {
                 mMultigridMap( Ik )( Ii ).set_size( tMaxIndex, 1, -1 );
            }
        }

        // Loop over all multigrid levels
        for ( moris::uint Ik = 0; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            for ( moris::uint Ij = 0; Ij < mListAdofExtIndMap( Ik ).length(); Ij++ )
            {
                moris::sint tTypeIdentifier = mListAdofTypeTimeIdentifier( Ik )( Ij, 0 );
                moris::sint tExtIndex = mListAdofExtIndMap( Ik )( Ij, 0 );

                mMultigridMap( Ik )( tTypeIdentifier )( tExtIndex, 0 ) = Ij;
            }
        }

        //////////////////////// printing/////////////////////
//        for ( moris::uint Ik = 0; Ik < mListAdofExtIndMap.size(); Ik++ )
//        {
//            print( mListAdofExtIndMap( Ik ), "mListAdofExtIndMap Ik");
//
//            for ( moris::uint Ij = 0; Ij < mMultigridMap( Ik ).size(); Ij++ )
//            {
//               print( mMultigridMap( Ik )( Ij ), "mMultigridMap Ik Ij");
//            }
//        }
    }

    //------------------------------------------------------------------------------------------------------------------------------------

    void Multigrid::read_multigrid_maps( const moris::uint               aLevel,
                                         const moris::Matrix< DDSMat > & aExtFineIndices,
                                         const moris::sint               aTypeTimeIdentifier,
                                               moris::Matrix< DDSMat > & aInternalFineIndices )
    {
        moris::uint tNumOfIndicesAsk = aExtFineIndices.length();

        aInternalFineIndices.set_size( tNumOfIndicesAsk, 1, -1 );

        for ( moris::uint Ik = 0; Ik < tNumOfIndicesAsk; Ik++ )
        {
            aInternalFineIndices( Ik, 0 ) = mMultigridMap( aLevel )( aTypeTimeIdentifier )( aExtFineIndices( Ik, 0 ), 0 );
        }
    }
}
}
