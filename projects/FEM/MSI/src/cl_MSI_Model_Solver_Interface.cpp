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

namespace moris
{
    namespace MSI
    {

        //------------------------------------------------------------------------------

        Model_Solver_Interface::Model_Solver_Interface(
                ParameterList                                       aMSIParameterList,
                std::shared_ptr< MSI::Equation_Model >              aEquationModel,
                mtk::Mesh                                         * aMesh )
        : mMSIParameterList( aMSIParameterList ),
          mEquationBlocks( aEquationModel->get_equation_sets() ),
          mDofMgn( aMesh->get_communication_table(), this ),
          mMesh( aMesh ),
          mEquationModel( aEquationModel )
        {
            this->create_equation_object_list();

            uint tNumCoeffs = 0;

            for( uint Ik = 0; Ik < mMesh->get_num_interpolations(); Ik++)
            {
                tNumCoeffs = std::max( tNumCoeffs, mMesh->get_max_num_coeffs_on_proc( Ik ) );
            }

            mDofMgn.set_max_num_adofs( tNumCoeffs );

            mDofMgn.initialize_pdof_type_list( mEquationBlocks );
        }

        //------------------------------------------------------------------------------

        void Model_Solver_Interface::finalize()
        {
            // get map from mesh
            Cell< moris::map< moris::moris_id, moris::moris_index > > & tCoefficientsIdtoIndexMap = mDofMgn.get_adof_map();

            if( mMesh != nullptr)         //FIXME fix all constructors to be able to get rid of this if statement
            {
                // get num discretizations
                uint tNumDiscretizations = mMesh->get_num_interpolations();

                // resize container for id to index maps
                tCoefficientsIdtoIndexMap.resize( tNumDiscretizations );

                // let mesh fill id to index maps
                for( uint Ik = 0; Ik < tNumDiscretizations; Ik++ )
                {
                    mMesh->get_adof_map( Ik, tCoefficientsIdtoIndexMap( Ik ) );
                }
            }

            for ( luint Ik = 0; Ik < mEquationBlocks.size(); ++Ik )
            {
                mEquationBlocks( Ik )->set_model_solver_interface( this );
            }

            mDofMgn.initialize_pdof_host_list( mEquationObjectList );

            mDofMgn.create_adofs();

            mDofMgn.set_pdof_t_matrix();

            for ( Equation_Object* tElement : mEquationObjectList )
            {
                tElement->create_my_pdof_list();

                tElement->create_my_list_of_adof_ids();

                tElement->set_unique_adof_map();
            }

            if ( mMSIParameterList.get< bool >( "multigrid" ) )
            {
                mMultigrid = new Multigrid( this, mMesh );

                mMultigrid->multigrid_initialize();
            }

            if( mMSIParameterList.get< bool >( "msi_checker" ) )
            {
                this->msi_checker();
            }
        }

        //------------------------------------------------------------------------------

        moris::sint Model_Solver_Interface::get_max_adof_index()
        {
            // max BSpline mesh index. However, not limited to BSpline meshes.
            sint tMaxAdofIndex = 0;

            uint tNumDofTypes = mDofMgn.get_num_dof_types();

            // loop over all used Ddf type indices
            for( uint Ik = 0; Ik< tNumDofTypes; Ik++)
            {
                sint tIndex = this->get_adof_index_for_type( Ik );

                tMaxAdofIndex = std::max( tIndex, tMaxAdofIndex );
            }

            return tMaxAdofIndex+1;
        }

        //------------------------------------------------------------------------------

        moris::sint Model_Solver_Interface::get_adof_index_for_type( moris::uint aDofType )
        {
            // Note: Make sure to add for each DOF type a default interpolation index to the MSI Parameter list

            // Get dof type enum
            enum Dof_Type tDofType = mDofMgn.get_dof_type_enum( aDofType );

            // get the index for the DoF-Type enum
            return this->get_adof_index_for_type( tDofType );
        }

        //------------------------------------------------------------------------------

        moris_index Model_Solver_Interface::get_adof_index_for_type( MSI::Dof_Type aDofType )
        {
            // Note: Make sure to add for each DOF type a default interpolation index to the MSI Parameter list

            if      ( aDofType == Dof_Type::TEMP )        { return mMSIParameterList.get< moris::sint >( "TEMP" ); }
            else if ( aDofType == Dof_Type::P )           { return mMSIParameterList.get< moris::sint >( "P" );    }
            else if ( aDofType == Dof_Type::RHO )         { return mMSIParameterList.get< moris::sint >( "RHO" );  }
            else if ( aDofType == Dof_Type::E )           { return mMSIParameterList.get< moris::sint >( "E" );    }
            else if ( aDofType == Dof_Type::EVP )         { return mMSIParameterList.get< moris::sint >( "EVP" );  }
            else if ( aDofType == Dof_Type::EVT )         { return mMSIParameterList.get< moris::sint >( "EVT" );  }

            else if ( aDofType == Dof_Type::UX )          { return mMSIParameterList.get< moris::sint >( "UX" );   }
            else if ( aDofType == Dof_Type::UY )          { return mMSIParameterList.get< moris::sint >( "UY" );   }
            else if ( aDofType == Dof_Type::UZ )          { return mMSIParameterList.get< moris::sint >( "UZ" );   }
            else if ( aDofType == Dof_Type::VX )          { return mMSIParameterList.get< moris::sint >( "VX" );   }
            else if ( aDofType == Dof_Type::VY )          { return mMSIParameterList.get< moris::sint >( "VY" );   }
            else if ( aDofType == Dof_Type::VZ )          { return mMSIParameterList.get< moris::sint >( "VZ" );   }
            else if ( aDofType == Dof_Type::MX )          { return mMSIParameterList.get< moris::sint >( "MX" );   }
            else if ( aDofType == Dof_Type::MY )          { return mMSIParameterList.get< moris::sint >( "MY" );   }
            else if ( aDofType == Dof_Type::MZ )          { return mMSIParameterList.get< moris::sint >( "MZ" );   }
            else if ( aDofType == Dof_Type::EVX )         { return mMSIParameterList.get< moris::sint >( "EVX" );  }
            else if ( aDofType == Dof_Type::EVY )         { return mMSIParameterList.get< moris::sint >( "EVY" );  }
            else if ( aDofType == Dof_Type::EVZ )         { return mMSIParameterList.get< moris::sint >( "EVZ" );  }
            else if ( aDofType == Dof_Type::NLSX )        { return mMSIParameterList.get< moris::sint >( "NLSX" ); }
            else if ( aDofType == Dof_Type::NLSY )        { return mMSIParameterList.get< moris::sint >( "NLSY" ); }
            else if ( aDofType == Dof_Type::NLSZ )        { return mMSIParameterList.get< moris::sint >( "NLSZ" ); }

            else if ( aDofType == Dof_Type::L2 )          { return mMSIParameterList.get< moris::sint >( "L2" );  }
            else if ( aDofType == Dof_Type::MAPPING_DOF ) { return mMSIParameterList.get< moris::sint >( "MAPPING_DOF" ); }
            else if ( aDofType == Dof_Type::LS1 )         { return mMSIParameterList.get< moris::sint >( "LS1" ); }
            else if ( aDofType == Dof_Type::LS2 )         { return mMSIParameterList.get< moris::sint >( "LS2" ); }

            else if ( aDofType == Dof_Type::THETA )       { return mMSIParameterList.get< moris::sint >( "THETA" ); }
            else if ( aDofType == Dof_Type::PHID )        { return mMSIParameterList.get< moris::sint >( "PHID" );  }
            else if ( aDofType == Dof_Type::PHISD )       { return mMSIParameterList.get< moris::sint >( "PHISD" ); }

            else if ( aDofType == Dof_Type::VISCOSITY )   { return mMSIParameterList.get< moris::sint >( "VISCOSITY" );  }
            else if ( aDofType == Dof_Type::STRESS_DOF )  { return mMSIParameterList.get< moris::sint >( "STRESS_DOF" ); }

            else if ( aDofType == Dof_Type::UNDEFINED )   { return MORIS_INDEX_MAX; }

            else
            {
                MORIS_ERROR( false, "Model_Solver_Interface::get_adof_index_for_type(): Dof type does not exist. Check dof type enums");
                return 0;
            }
        }

        //------------------------------------------------------------------------------

        void Model_Solver_Interface::msi_checker()
        {
            this->pdof_host_checker();

            // add more msi checks here!!!
        }

        void Model_Solver_Interface::pdof_host_checker()
        {
            for ( luint Ik = 0; Ik < mEquationBlocks.size(); ++Ik )
            {
                Cell< MSI::Equation_Object * > & tEquationObj = mEquationBlocks( Ik )->get_equation_object_list();

                // variables used for number of system check.
                uint tNumEquationSys = 0;
                bool tNumEquationSysSet = false;

                // variables used for pdof host check
                Matrix< DDUMat > tNumPdofsOnPdofHosts;
                uint tNumPotentialPdofs = 0;
                bool tNumPdofsOnPdofHostsSet = false;

                for ( luint Ij = 0; Ij < tEquationObj.size(); ++Ij )
                {
                    const moris::Cell< moris::Cell< Pdof_Host * > > & tPdofHosts = tEquationObj( Ij )->get_pdof_hosts();

                    if( !tNumEquationSysSet )
                    {
                        tNumEquationSysSet = true;
                        tNumEquationSys = tPdofHosts.size();
                    }

                    MORIS_ERROR( tPdofHosts.size() == tNumEquationSys,
                            "Equation Objects on Set do not have the same number of systems. This might indicate that there are elements without vertices");

                    for ( luint Ii = 0; Ii < tPdofHosts.size(); ++Ii )
                    {
                        for ( luint Ia = 0; Ia < tPdofHosts( Ii ).size(); ++Ia )
                        {
                            moris::Cell< moris::Cell< Pdof* > > & tPdofHostPdofList = tPdofHosts( Ii )( Ia )-> get_pdof_hosts_pdof_list();

                            if( !tNumPdofsOnPdofHostsSet )
                            {
                                tNumPdofsOnPdofHostsSet = true;
                                tNumPotentialPdofs = tPdofHostPdofList.size();
                                tNumPdofsOnPdofHosts.set_size( tNumPotentialPdofs,1, 0 );

                                for ( luint Ib = 0; Ib < tNumPotentialPdofs; ++Ib )
                                {
                                    tNumPdofsOnPdofHosts( Ib ) = tPdofHostPdofList( Ib ).size();
                                }
                            }

                            MORIS_ERROR( tNumPotentialPdofs == tPdofHostPdofList.size(),
                                    "Some Pdof hosts on set have different number of potential pdofs. Check global dof type list!");

                            for ( luint Ib = 0; Ib < tPdofHostPdofList.size(); ++Ib )
                            {
                                MORIS_ERROR( tNumPdofsOnPdofHosts( Ib ) == tPdofHostPdofList( Ib ).size(),
                                        "Some Pdof hosts on set do not share the same dofs. Check dof dependencies in input file. (Especially ghost in parallel)");
                            }
                        }
                    }
                }
            }
        }
    }
}

