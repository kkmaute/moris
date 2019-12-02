/*
 * cl_FEM_IWG.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: sonne
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"

namespace moris
{
namespace fem
{

void IWG::get_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
{
    // set the size of the dof type list for the set
    uint tCounter = 0;

    for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
    {
        tCounter += mMasterDofTypes( iDOF ).size();
    }
    for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
    {
        tCounter += mSlaveDofTypes( iDOF ).size();
    }

    for ( std::shared_ptr< Property > tProperty : mMasterProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    for ( std::shared_ptr< Property > tProperty : mSlaveProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type >  tActiveDofType;
        tCM->get_dof_types( tActiveDofType );

        tCounter += tActiveDofType.size();
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

        tCounter += tActiveDofType.size();
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

        tCounter += tActiveDofType.size();
    }

    aDofTypes.reserve( tCounter );

    for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
    {
        aDofTypes.append( mMasterDofTypes( iDOF ) );
    }
    for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
    {
        aDofTypes.append( mSlaveDofTypes( iDOF )  );
    }

    for ( std::shared_ptr< Property > tProperty : mMasterProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    for ( std::shared_ptr< Property > tProperty : mSlaveProp )
    {
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tProperty->get_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from constitutive models
    for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
    {
        // get dof types for constitutive model
        moris::Cell< MSI::Dof_Type > tActiveDofType;
        tCM->get_dof_types( tActiveDofType );

        aDofTypes.append( tActiveDofType );
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

        // loop on property dof type
        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
        {
            aDofTypes.append( tActiveDofType( iDOF ) );
        }
    }

    // get dof type from stabilization parameters
    for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
    {
        // get dof types for constitutive model
        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

        // loop on property dof type
        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
        {
            aDofTypes.append( tActiveDofType( iDOF ) );
        }
    }
}

void IWG::build_global_dof_type_list()
    {
        // MASTER-------------------------------------------------------
        // set size for the global dof type list
        uint tNumDofTypes = mSet->get_num_dof_types();

        // set size for the global dof type list
        mMasterGlobalDofTypes.reserve( tNumDofTypes );

        // set a size for the checkList (used to avoid repeating a dof type)
        Matrix< DDSMat > tCheckList( tNumDofTypes, 1, -1 );

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
        {
            sint tDofTypeindex = mSet->get_dof_index_for_type( mMasterDofTypes( iDOF )( 0 ) );  //FIXME'

            // put the dof type in the checklist
            tCheckList( tDofTypeindex ) = 1;

            // put the dof type in the global type list
            mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDOF ) );
        }

        // get dof type from properties
        for ( std::shared_ptr< Property > tProperty : mMasterProp )
        {
            // get dof types for property
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );  //FIXME

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }

        // get dof type from constitutive models
        for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }

        // get dof type from stabilization parameters
        for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }

        // SLAVE--------------------------------------------------------

        // set size for the global dof type list
        mSlaveGlobalDofTypes.reserve( tNumDofTypes );

        // set a size for the checkList (used to avoid repeating a dof type)
        tCheckList.fill( -1 );

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
        {
            sint tDofTypeindex = mSet->get_dof_index_for_type( mSlaveDofTypes( iDOF )( 0 ) );

            // put the dof type in the checklist
            tCheckList( tDofTypeindex ) = 1;

            // put the dof type in the global type list
            mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDOF ) );
        }

        // get dof type from properties
        for ( std::shared_ptr< Property > tProperty : mSlaveProp )
        {
            // get dof types for property
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }

        // get dof type from constitutive models
        for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }

        // get dof type from stabilization parameters
        for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
        {
            // get dof types for constitutive model
            moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

            // loop on property dof type
            for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
            {
                sint tDofTypeindex = mSet->get_dof_index_for_type( tActiveDofType( iDOF )( 0 ) );

                // if dof enum not in the list
                if ( tCheckList( tDofTypeindex) != 1 )
                {
                    // put the dof type in the checklist
                    tCheckList( tDofTypeindex ) = 1;

                    // put the dof type in the global type list
                    mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                }
            }
        }
    }
//------------------------------------------------------------------------------

}   // end fem namespace
}   // end moris namespace


