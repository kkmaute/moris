/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Model_Solver_Interface.cpp
 *
 */

#include "cl_MSI_Model_Solver_Interface.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_MSI_get_mesh_index_for_dof_type.hpp"

namespace moris::MSI
{

    //------------------------------------------------------------------------------

    Model_Solver_Interface::Model_Solver_Interface(
            const Parameter_List&                         aMSIParameterList,
            const std::shared_ptr< MSI::Equation_Model >& aEquationModel,
            mtk::Mesh*                                    aMesh )
            : mMSIParameterList( aMSIParameterList )
            , mEquationSets( aEquationModel->get_equation_sets() )
            , mDofMgn( aMesh->get_communication_table(), this )
            , mMesh( aMesh )
            , mEquationModel( aEquationModel )
    {
        this->create_equation_object_list();

        // determine maximum number of generic adof among all interpolation mesh
        uint tNumCoeffs = 0;

        for ( uint Ik = 0; Ik < mMesh->get_num_interpolations(); Ik++ )
        {
            tNumCoeffs = std::max( tNumCoeffs, mMesh->get_max_num_coeffs_on_proc( Ik ) );
        }

        // set maximum number of generic adof among all interpolation meshes
        mDofMgn.set_max_num_adofs( tNumCoeffs );

        // create list of all pdof types used in equation objects
        mDofMgn.initialize_pdof_type_list( mEquationSets );
    }

    //------------------------------------------------------------------------------

    void
    Model_Solver_Interface::create_equation_object_list()
    {
        // clear list but do not delete equation object pointer
        mEquationObjectList.clear();

        // reset number of equation object count
        moris::uint tNumEquationObj = 0;

        for ( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
        {
            tNumEquationObj += mEquationSets( Ik )->get_num_equation_objects();
        }

        mEquationObjectList.reserve( tNumEquationObj );

        for ( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
        {
            mEquationObjectList.append( mEquationSets( Ik )->get_equation_object_list() );
        }
    }

    //------------------------------------------------------------------------------

    void
    Model_Solver_Interface::finalize()
    {
        // get map from mesh
        Vector< moris::map< moris::moris_id, moris::moris_index > >& tCoefficientsIdtoIndexMap = mDofMgn.get_adof_map();

        if ( mMesh != nullptr )    // FIXME fix all constructors to be able to get rid of this if statement
        {
            // get num discretizations
            uint tNumDiscretizations = mMesh->get_num_interpolations();

            // resize container for id to index maps
            tCoefficientsIdtoIndexMap.resize( tNumDiscretizations );

            // let mesh fill id to index maps
            for ( uint Ik = 0; Ik < tNumDiscretizations; Ik++ )
            {
                mMesh->get_adof_map( Ik, tCoefficientsIdtoIndexMap( Ik ) );
            }
        }

        for ( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
        {
            mEquationSets( Ik )->set_model_solver_interface( this );
        }

        mDofMgn.initialize_pdof_host_list( mEquationObjectList );

        // create adofs
        mDofMgn.create_adofs();

        // set T matrix
        mDofMgn.set_pdof_t_matrix();

        for ( Equation_Object* tEquationObject : mEquationObjectList )
        {
            tEquationObject->create_my_pdof_list();

            tEquationObject->create_my_list_of_adof_ids();

            tEquationObject->set_unique_adof_map();
        }

        if ( mMSIParameterList.get< bool >( "multigrid" ) )
        {
            mMultigrid = new Multigrid( this, mMesh );

            mMultigrid->multigrid_initialize();
        }

        if ( mMSIParameterList.get< bool >( "msi_checker" ) )
        {
            this->msi_checker();
        }
    }

    //------------------------------------------------------------------------------

    moris::sint
    Model_Solver_Interface::get_max_adof_index()
    {
        // max BSpline mesh index. However, not limited to BSpline meshes.
        sint tMaxAdofIndex = 0;

        uint tNumDofTypes = mDofMgn.get_num_dof_types();

        // loop over all used Ddf type indices
        for ( uint Ik = 0; Ik < tNumDofTypes; Ik++ )
        {
            sint tIndex = this->get_adof_index_for_type( Ik );

            tMaxAdofIndex = std::max( tIndex, tMaxAdofIndex );
        }

        return tMaxAdofIndex + 1;
    }

    //------------------------------------------------------------------------------

    moris::sint
    Model_Solver_Interface::get_adof_index_for_type( moris::uint aDofType )
    {
        // Note: Make sure to add for each DOF type a default interpolation index to the MSI Parameter list

        // Get dof type enum
        enum Dof_Type tDofType = mDofMgn.get_dof_type_enum( aDofType );

        // get the index for the DoF-Type enum
        return this->get_adof_index_for_type( tDofType );
    }

    //------------------------------------------------------------------------------

    moris_index
    Model_Solver_Interface::get_adof_index_for_type( MSI::Dof_Type aDofType )
    {
        // just call the designated function and pass the MSI parameter list
        return get_mesh_index_for_dof_type( aDofType, mMSIParameterList );
    }

    //------------------------------------------------------------------------------

    void
    Model_Solver_Interface::msi_checker()
    {
        this->pdof_host_checker();

        // add more msi checks here!!!
    }

    //------------------------------------------------------------------------------

    void
    Model_Solver_Interface::pdof_host_checker()
    {
        for ( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
        {
            Vector< MSI::Equation_Object* >& tEquationObj = mEquationSets( Ik )->get_equation_object_list();

            // variables used for number of system check.
            uint tNumEquationSys    = 0;
            bool tNumEquationSysSet = false;

            // variables used for pdof host check
            Matrix< DDUMat > tNumPdofsOnPdofHosts;
            uint             tNumPotentialPdofs      = 0;
            bool             tNumPdofsOnPdofHostsSet = false;

            for ( luint Ij = 0; Ij < tEquationObj.size(); ++Ij )
            {
                const Vector< Vector< Pdof_Host* > >& tPdofHosts = tEquationObj( Ij )->get_pdof_hosts();

                if ( !tNumEquationSysSet )
                {
                    tNumEquationSysSet = true;
                    tNumEquationSys    = tPdofHosts.size();
                }

                MORIS_ERROR( tPdofHosts.size() == tNumEquationSys,
                        "Equation Objects on Set do not have the same number of systems. This might indicate that there are elements without vertices" );

                for ( luint Ii = 0; Ii < tPdofHosts.size(); ++Ii )
                {
                    for ( luint Ia = 0; Ia < tPdofHosts( Ii ).size(); ++Ia )
                    {
                        Vector< Vector< Pdof* > >& tPdofHostPdofList = tPdofHosts( Ii )( Ia )->get_pdof_hosts_pdof_list();

                        if ( !tNumPdofsOnPdofHostsSet )
                        {
                            tNumPdofsOnPdofHostsSet = true;
                            tNumPotentialPdofs      = tPdofHostPdofList.size();
                            tNumPdofsOnPdofHosts.set_size( tNumPotentialPdofs, 1, 0 );

                            for ( luint Ib = 0; Ib < tNumPotentialPdofs; ++Ib )
                            {
                                tNumPdofsOnPdofHosts( Ib ) = tPdofHostPdofList( Ib ).size();
                            }
                        }

                        MORIS_ERROR( tNumPotentialPdofs == tPdofHostPdofList.size(),
                                "Some Pdof hosts on set have different number of potential pdofs. Check global dof type list!" );

                        for ( luint Ib = 0; Ib < tPdofHostPdofList.size(); ++Ib )
                        {
                            MORIS_ERROR( tNumPdofsOnPdofHosts( Ib ) == tPdofHostPdofList( Ib ).size(),
                                    "Some Pdof hosts on set do not share the same dofs. Check dof dependencies in input file. (Especially ghost in parallel)" );
                        }
                    }
                }
            }
        }
    }
}    // namespace moris::MSI
