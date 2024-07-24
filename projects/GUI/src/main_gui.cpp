/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main_gui.cpp
 *
 */


#include "main_gui.hpp"


// moris::Comm_Manager gMorisComm;
// moris::Logger       gLogger;
namespace moris
{
    Moris_Gui::Moris_Gui( QWidget *parent )
            : QWidget( parent )
    {
        mParameterLists = read();
        // mLayout is the main layout of the GUI
        // mSidePanel is the side panel layout where the tree widget and buttons are placed
        mLayout->addLayout( mSidePanel );
        mTreeWidget->setColumnCount( 1 );

        // Read from mParamterLists instead of moris functions, put these functions in the read() function

        for ( uint iProject = 0; iProject < (uint)( moris::Parameter_List_Type::END_ENUM ); iProject++ )
        {
            mProjectNames.append( QString::fromStdString( moris::convert_parameter_list_enum_to_string( ( moris::Parameter_List_Type )( iProject ) ) ) );
        }

        moris::Parameter_List_Type tModule = moris::Parameter_List_Type::OPT;
        QStringList                tStringList;
        QStringList                tSubChildrenList;

        // Initialize data structures for each project
        mTreeWidgetItems.resize( mProjectNames.size() );
        mTreeWidgetChildren.resize( mProjectNames.size() );
        mTreeWidgetSubChildren.resize( mProjectNames.size() );

        // QTreeWidgetItems
        mQTreeWidgetItems.resize( mProjectNames.size() );
        mQTreeWidgetChildren.resize( mProjectNames.size() );
        mQTreeWidgetSubChildren.resize( mProjectNames.size() );


        for ( uint iRoot = 0; iRoot < mProjectNames.size(); iRoot++ )
        {
            // Get the number of sub-parameter lists in the module
            tStringList.clear();

            // Read from mParamterLists instead of moris functions, put these functions in the read() function
            for ( uint iList = 0; iList < moris::get_number_of_sub_parameter_lists_in_module( tModule ); iList++ )
            {
                tStringList.append( QString::fromStdString( moris::get_outer_sub_parameter_list_name( tModule, iList ) ) );
            }

            // Create the tree widget items for the project names
            mQTreeWidgetItems[ iRoot ] = new QTreeWidgetItem( mTreeWidget );
            mQTreeWidgetItems[ iRoot ]->setText( 0, mProjectNames[ iRoot ] );

            mTreeWidgetItems[ iRoot ] = new Moris_Tree_Widget_Item();

            // Resize the children and sub-children lists
            mTreeWidgetChildren[ iRoot ].resize( tStringList.size() );
            mTreeWidgetSubChildren[ iRoot ].resize( tStringList.size() );

            mQTreeWidgetChildren[ iRoot ].resize( tStringList.size() );
            mQTreeWidgetSubChildren[ iRoot ].resize( tStringList.size() );

            // For each of the children of the project,
            for ( uint iChildren = 0; iChildren < tStringList.size(); iChildren++ )
            {

                // Set the TreeWidgetItems for the children and set their texts
                mQTreeWidgetChildren[ iRoot ][ iChildren ] = new QTreeWidgetItem( mQTreeWidgetItems[ iRoot ] );
                mQTreeWidgetChildren[ iRoot ][ iChildren ]->setText( 0, tStringList[ iChildren ] );

                mTreeWidgetChildren[ iRoot ][ iChildren ] = new Moris_Tree_Widget_Item();

                // Set up the scroll area for the form layouts as these will have associated forms
                mTreeWidgetChildren[ iRoot ][ iChildren ]->setupScrollArea();

                // Storing the index of the Moris_Tree_Widget_Item
                mTreeWidgetChildren[ iRoot ][ iChildren ]->setIndex( { iRoot, iChildren } );

                if ( iRoot == (uint)( moris::Parameter_List_Type::OPT ) && iChildren == (uint)( moris::OPT_SubModule::OPTIMIZATION_PROBLEMS ) )
                {
                    // Setting the OPT/Optimization Problems form visible upon construction of the GUI
                    std::cout << mParameterLists( iRoot )( iChildren ).size() << std::endl;

                    mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mParameterLists( iRoot )( iChildren )( 0 ) );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( true );
                    mOldItem = mQTreeWidgetChildren[ iRoot ][ iChildren ];
                }
                else if ( iRoot == (uint)( moris::Parameter_List_Type::OPT ) && iChildren == (uint)( moris::OPT_SubModule::ALGORITHMS ) )
                {
                    // Adding the OPT/Algorithms to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( { "gcmma", "lbfgs", "sql", "sweep" } );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );
                }
                else if ( iRoot == (uint)( moris::Parameter_List_Type::GEN ) && iChildren == (uint)( moris::GEN_SubModule::GEOMETRIES ) )
                {
                    // Adding the GEN/Geometries to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    QStringList tGENGeometriesList;
                    for ( uint iGeometries = 0; iGeometries < (uint)moris::gen::Geometry_Type_String::values.size(); iGeometries++ )
                    {
                        // For Level Set geometries, iterates through the Field_Types
                        if ( (uint)iGeometries == (uint)moris::gen::Geometry_Type::LEVEL_SET )
                        {
                            for ( uint iLevelSet = 0; iLevelSet < (uint)moris::gen::Field_Type_String::values.size(); iLevelSet++ )
                            {
                                tGENGeometriesList.append( QString::fromStdString( moris::gen::Field_Type_String::values( iLevelSet ) ) );
                            }
                        }
                        else
                        {
                            // For other geometries, iterates through the Geometry_Types
                            tGENGeometriesList.append( QString::fromStdString( moris::gen::Geometry_Type_String::values( iGeometries ) ) );
                        }
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tGENGeometriesList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );
                }
                else if ( iRoot == (uint)( moris::Parameter_List_Type::SOL ) && iChildren == (uint)( moris::SOL_SubModule::LINEAR_ALGORITHMS ) )
                {
                    // Adding the SOL/Linear Algorithms to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    QStringList tSolverList;
                    for ( uint iSolver = 0; iSolver < (uint)moris::sol::SolverType_String::values.size() - 1; iSolver++ )
                    {
                        tSolverList.append( QString::fromStdString( moris::sol::SolverType_String::values( iSolver ) ) );
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tSolverList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );
                }
                else
                {
                    // Adding the parameter_list elements to all the forms and setting the visibility to false
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mParameterLists( iRoot )( iChildren )( 0 ) );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );
                }

                // Adding the scroll area to the main layout
                mLayout->addWidget( mTreeWidgetChildren[ iRoot ][ iChildren ]->getScrollArea() );

                // Adds the children to the associated project
                mQTreeWidgetItems[ iRoot ]->addChild( mQTreeWidgetChildren[ iRoot ][ iChildren ] );

                // Connects the QTreeWidgetChildren to the associated Moris_Tree_Widget_Item
                mTreeWidget->setItemWidget( mQTreeWidgetChildren[ iRoot ][ iChildren ], 0, mTreeWidgetChildren[ iRoot ][ iChildren ] );
            }

            // Adds the root to the tree widget
            mTreeWidget->addTopLevelItem( mQTreeWidgetItems[ iRoot ] );

            // Connects the QTreeWidgetItems to the associated Moris_Tree_Widget_Item
            mTreeWidget->setItemWidget( mQTreeWidgetItems[ iRoot ], 0, mTreeWidgetItems[ iRoot ] );

            tModule = ( moris::Parameter_List_Type )( (uint)( tModule ) + 1 );

            // MORIS_ERROR( false, "error here" );
        }

        // Setting the header label for the tree widget
        mTreeWidget->setHeaderLabel( "Projects" );

        // Adding the tree widget to the side panel
        mSidePanel->addWidget( mTreeWidget );

        // Setting the text for the add and remove buttons
        mAddButton->setText( "Add" );
        mSidePanel->addWidget( mAddButton );

        mRemoveButton->setText( "Remove" );
        mSidePanel->addWidget( mRemoveButton );

        // Connecting the signals and slots
        connect( mTreeWidget, SIGNAL( itemSelectionChanged() ), this, SLOT( parameter_selected() ) );
        connect( mAddButton, SIGNAL( clicked() ), this, SLOT( add_more_props() ) );
        connect( mRemoveButton, SIGNAL( clicked() ), this, SLOT( remove_props() ) );
    }

    void Moris_Gui::parameter_selected()
    {

        /**
         * @brief Function to keep track of the selected item in the tree widget
         * @param NONE
         * @return NONE
         * @note This function keeps track of the selected item in the tree widget and sets the form associated with the selected item to visible.
         */

        QTreeWidgetItem *tCurrentItem = mTreeWidget->currentItem();
        if ( tCurrentItem )
        {
            Moris_Tree_Widget_Item *tCurrentWidget = qobject_cast< Moris_Tree_Widget_Item * >( mTreeWidget->itemWidget( tCurrentItem, 0 ) );
            if ( tCurrentWidget->hasForm() )
            {
                Moris_Tree_Widget_Item *tOldWidget = qobject_cast< Moris_Tree_Widget_Item * >( mTreeWidget->itemWidget( mOldItem, 0 ) );
                tOldWidget->set_form_visible( false );
                tCurrentWidget->set_form_visible( true );
                mOldItem = tCurrentItem;
            }
        }
    }


    QList< QStringList > Moris_Gui::convert_parameters_to_QStringList( moris::Parameter_List aList )
    {

        /**
         * @brief Function to convert the parameter list to a QList of QStringList
         * @param moris::Parameter_List aParameterList
         * @return QList< QStringList >
         * @note This function converts the MORIS parameter list to a QList of QStringList where the 0th index of the list gives the key of the parameter_list and the 1st index gives the value of the parameter_list
         * @note tQtList[0] corresponds to the "key" values of the parameter map, tQtList[1] corresponds to the "value" of the parameter map
         */

        QList< QStringList > tQtList;
        QStringList          tStringList;
        tQtList.append( tStringList );
        tQtList.append( tStringList );

        for ( auto it = aList.begin(); it != aList.end(); ++it )
        {
            tQtList[ 0 ].append( QString::fromStdString( it->first ) );
            tQtList[ 1 ].append( QString::fromStdString( it->second.get_string() ) );
        }

        return tQtList;
    }

    void Moris_Gui::add_more_props()
    {
        /**
         * @brief Function to add more parameters when the mAddButton is clicked
         * @note This function is called when the mAddButton is clicked. It adds more parameters in a form or adds a new sub-form if the Moris_Tree_Widget_Item has a sub-form.
         */

        QTreeWidgetItem *tCurrentItem = mTreeWidget->currentItem();

        // Check if the current item is not null
        if ( tCurrentItem )
        {

            // Making sure there are no casting problems when getting the current widget
            Moris_Tree_Widget_Item *tCurrentWidget = qobject_cast< Moris_Tree_Widget_Item * >( mTreeWidget->itemWidget( tCurrentItem, 0 ) );

            // Check if the current widget has a sub-form
            if ( tCurrentWidget->hasSubForm() )
            {

                // Get the index of the current widget
                QList< uint > tCurrentIndex = tCurrentWidget->getIndex();

                // Append a new QTreeWidgetItem and Moris_Tree_Widget_Item to the sub-children list of the current widget
                // And set up the necessary elements
                qDebug() << "Adding Sub Form";
                mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].append( new QTreeWidgetItem( mQTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ] ) );
                mQTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->addChild( mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last() );
                mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setText( 0, tCurrentWidget->getComboBox()->currentText() + QString::number( tCurrentWidget->getSubFormCountProps() ) );

                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].append( new Moris_Tree_Widget_Item() );

                // Reorganizing the count of the sub-forms from the item selection in the combo box
                mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->addSubFormCountProps();

                // Setting up the scroll area, index, saving type for the sub-forms
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setupScrollArea();
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setIndex( { tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], (uint)mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].size() - 1 } );
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setSubFormType( tCurrentWidget->getComboBox()->currentIndex() );
                // Append to mParameterLists and pass it into add_elements
                mParameterLists( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).push_back( create_parameter_list( (moris::Parameter_List_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], tCurrentWidget->getComboBox()->currentIndex() ) );
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->add_elements( mParameterLists( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).back() );
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->set_form_visible( false );

                // Connecting the QTreeWidgetItem to the Moris_Tree_Widget_Item
                mTreeWidget->setItemWidget( mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last(), 0, mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last() );

                // Adding the scroll area to the main layout
                mLayout->addWidget( mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->getScrollArea() );
            }
            else
            {
                QList< uint > tCurrentIndex = tCurrentWidget->getIndex();
                tCurrentWidget->add_elements( mParameterLists( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] )(0) );
            }
        }
    }

    void Moris_Gui::remove_props()
    {
        /**
         * @brief Function to remove parameters when the mRemoveButton is clicked
         * @note This function is called when the mRemoveButton is clicked. It removes parameters in a form or removes a sub-form if the Moris_Tree_Widget_Item has sub-forms.
         */

        QTreeWidgetItem *tCurrentItem = mTreeWidget->currentItem();

        // Check if the current item is not null
        if ( tCurrentItem )
        {

            // Making sure there are no casting problems when getting the current widget
            Moris_Tree_Widget_Item *tCurrentWidget = qobject_cast< Moris_Tree_Widget_Item * >( mTreeWidget->itemWidget( tCurrentItem, 0 ) );

            // Check if the current widget has a sub-form
            if ( tCurrentWidget->isSubForm() )
            {
                // Get the index of the current widget
                QList< uint > tCurrentIndex = tCurrentWidget->getIndex();

                // Reorganizing the count of the sub-forms from the item selection in the combo box
                mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->removeSubFormCountProps( tCurrentWidget->getSubFormType() );

                // If the first form is selected then the user is redirected to the main form
                if ( tCurrentIndex[ 2 ] == 0 )
                {
                    mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->set_form_visible( true );
                }
                else
                {
                    mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ][ tCurrentIndex[ 2 ] ]->set_form_visible( true );
                }

                // Removing the sub-form from the layout and deleting the associated objects
                delete mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ][ tCurrentIndex[ 2 ] ];
                delete mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ][ tCurrentIndex[ 2 ] ];
                mParameterLists( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).erase( tCurrentIndex[ 2 ] );
            }
            else
            {
                // If the current widget is not a sub-form, then remove the elements from the form in the Moris_Tree_Widget_Item
                tCurrentWidget->remove_elements();
            }
        }
    }


    Vector< Vector< Vector< Parameter_List > > > read()
    {

        // Create the 3d vector
        Vector< Vector< Vector< Parameter_List > > > tParameterList;
        tParameterList.resize( (uint)( moris::Parameter_List_Type::END_ENUM ) );
        for ( uint iRoot = 0; iRoot < (uint)( moris::Parameter_List_Type::END_ENUM ); iRoot++ )
        {
            tParameterList( iRoot ).resize( moris::get_number_of_sub_parameter_lists_in_module( (moris::Parameter_List_Type)iRoot ) );
            for ( uint iChild = 0; iChild < moris::get_number_of_sub_parameter_lists_in_module( (moris::Parameter_List_Type)iRoot ); iChild++ )
            {
                if ( ( iRoot == (uint)( moris::Parameter_List_Type::OPT ) && iChild == (uint)( moris::OPT_SubModule::ALGORITHMS ) )
                        || ( iRoot == (uint)( moris::Parameter_List_Type::GEN ) && iChild == (uint)( moris::GEN_SubModule::GEOMETRIES ) )
                        || ( iRoot == (uint)( moris::Parameter_List_Type::SOL ) && iChild == (uint)( moris::SOL_SubModule::LINEAR_ALGORITHMS ) ) )
                {
                }
                else {
                    tParameterList( iRoot )( iChild ).push_back( create_parameter_list( (moris::Parameter_List_Type)iRoot, iChild, 0 ) );
                    std::cout << "Root: " << iRoot << " Child: " << iChild << " Size: " << tParameterList( iRoot )( iChild ).size() << std::endl;
                }
            }
        }

        // Resize based on the projects/sub-projects
        // Populate based on the parameter_list
        return tParameterList;
    }

    Parameter_List create_parameter_list( Parameter_List_Type aModule, uint aChild, uint aSubChild )
    {
        /*
        function name: get_parameter_list
        parameters:
          moris::Parameter_List_Type aModule (ENUM) -> this gives the project name
          uint aChild -> gives the child index
          uint aSubChild -> gives the Sub-Child (inner sub-module) index
        returns:
            QList <QStringList>
                the create_function returns a ParameterList object that is a type of map
                The 0th index of the QList gives the "keys" of the map
                The 1st index of the QList gives the default "values" of the map
        */


        moris::Parameter_List tParameterList;

        switch ( aModule )
        {
            case moris::Parameter_List_Type::OPT:
                switch ( aChild )
                {
                    case 0:
                        tParameterList = ( moris::prm::create_opt_problem_parameter_list() );
                        break;

                    case 1:
                        tParameterList = ( moris::prm::create_opt_interface_parameter_list() );

                        // Commented out the Interface manager for now

                        // switch ( aSubChild )
                        // {
                        //     case 0:
                        //         tParameterList = ( moris::prm::create_opt_interface_parameter_list() );

                        //         break;

                        //     case 1:
                        //         tParameterList = ( moris::prm::create_opt_interface_manager_parameter_list() );

                        //         break;

                        //     default:
                        //         break;
                        // }

                        break;

                    case 2:
                        switch ( aSubChild )
                        {
                            case 0:
                                tParameterList = ( moris::prm::create_gcmma_parameter_list() );

                                break;

                            case 1:
                                tParameterList = ( moris::prm::create_lbfgs_parameter_list() );
                                break;

                            case 2:
                                tParameterList = ( moris::prm::create_sqp_parameter_list() );

                                break;

                            case 3:
                                tParameterList = ( moris::prm::create_sweep_parameter_list() );

                                break;
                            default:
                                break;
                        }

                        break;

                    default:
                        break;
                }
                // Free
                break;

            case moris::Parameter_List_Type::HMR:
                tParameterList = ( moris::prm::create_hmr_parameter_list() );

                break;

            case moris::Parameter_List_Type::STK:
                tParameterList = ( moris::prm::create_stk_parameter_list() );
                break;

            case moris::Parameter_List_Type::XTK:
                tParameterList = ( moris::prm::create_xtk_parameter_list() );
                break;

            case moris::Parameter_List_Type::GEN:
                switch ( aChild )
                {
                    case 0:
                    {
                        tParameterList = ( moris::prm::create_gen_parameter_list() );
                        break;
                    }

                    case 1:
                    {
                        if ( aSubChild <= (uint)moris::gen::Field_Type::USER_DEFINED )
                        {
                            tParameterList = ( moris::prm::create_level_set_geometry_parameter_list( (moris::gen::Field_Type)aSubChild ) );
                        }
                        else if ( aSubChild == (uint)moris::gen::Field_Type::USER_DEFINED + 1 )
                        {
                            tParameterList = ( moris::prm::create_surface_mesh_geometry_parameter_list() );
                        }
                        else
                        {
                            tParameterList = ( moris::prm::create_voxel_geometry_parameter_list() );
                        }
                        break;
                    }
                    case 2:
                    {
                        tParameterList = ( moris::prm::create_gen_property_parameter_list( moris::gen::Field_Type::CONSTANT ) );
                        break;
                    }
                    default:
                    {
                        break;
                    }
                }

                break;

            case moris::Parameter_List_Type::FEM:
                /*
                 * Set of Dropdowns for tParameterList[0] (property_name in FEM)
                 * //Dropdown
                 * PropDensity, PropYoungs, PropPoisson,
                 * PropCTE, PropRefTemp, PropConductivity,
                 * PropCapacity, PropDirichlet, PropSelectX,
                 * PropSelectY, PropSelectZ, PropInnerPressureLoad,
                 * Pro#include <QApplication>
                 * pOuterPressureLoad, PropOuterTemperature
                 */
                switch ( aChild )
                {
                    case 0:
                        tParameterList = ( moris::prm::create_property_parameter_list() );
                        break;

                    case 1:
                        tParameterList = ( moris::prm::create_constitutive_model_parameter_list() );
                        break;

                    case 2:
                        tParameterList = ( moris::prm::create_stabilization_parameter_parameter_list() );
                        break;

                    case 3:
                        tParameterList = ( moris::prm::create_IWG_parameter_list() );
                        break;

                    case 4:
                        tParameterList = ( moris::prm::create_IQI_parameter_list() );
                        break;

                    case 5:
                        tParameterList = ( moris::prm::create_computation_parameter_list() );
                        break;

                    case 6:
                        tParameterList = ( moris::prm::create_fem_field_parameter_list() );
                        break;

                    case 7:
                        tParameterList = ( moris::prm::create_phase_parameter_list() );
                        break;

                    case 8:
                        tParameterList = ( moris::prm::create_material_model_parameter_list() );
                        break;

                    default:
                        break;
                }


                break;

            case moris::Parameter_List_Type::SOL:
                //            tParameterList.resize( 8 );

                switch ( aChild )
                {
                    case 0:

                        switch ( aSubChild )
                        {
                            case 0:
                                tParameterList = ( moris::prm::create_linear_algorithm_parameter_list_aztec() );
                                break;

                            case 1:
                                tParameterList = ( moris::prm::create_linear_algorithm_parameter_list_amesos() );
                                break;

                            case 2:
                                tParameterList = ( moris::prm::create_linear_algorithm_parameter_list_belos() );
                                break;

                            case 3:
                                tParameterList = ( moris::prm::create_linear_algorithm_parameter_list_petsc() );
                                break;

                            case 4:
                                tParameterList = ( moris::prm::create_eigen_algorithm_parameter_list() );
                                break;

                            case 5:
                                // Need to add ML here
                                tParameterList = ( moris::prm::create_linear_algorithm_parameter_list_belos() );
                                break;

                            case 6:
                                tParameterList = ( moris::prm::create_slepc_algorithm_parameter_list() );
                                break;

                            default:
                                break;
                        }
                        break;

                        break;

                    case 1:
                        tParameterList = ( moris::prm::create_linear_solver_parameter_list() );
                        break;

                    case 2:
                        tParameterList = ( moris::prm::create_nonlinear_algorithm_parameter_list() );
                        break;

                    case 3:
                        tParameterList = ( moris::prm::create_nonlinear_solver_parameter_list() );
                        break;

                    case 4:
                        tParameterList = ( moris::prm::create_time_solver_algorithm_parameter_list() );
                        break;

                    case 5:
                        tParameterList = ( moris::prm::create_time_solver_parameter_list() );
                        break;

                    case 6:
                        tParameterList = ( moris::prm::create_solver_warehouse_parameterlist() );
                        break;

                    case 7:
                        // Need to add Preconditioners
                        tParameterList = ( moris::prm::create_material_model_parameter_list() );
                        break;

                    default:
                        break;
                }

                break;

            case moris::Parameter_List_Type::MSI:
                tParameterList = ( moris::prm::create_msi_parameter_list() );
                break;

            case moris::Parameter_List_Type::VIS:
                tParameterList = ( moris::prm::create_vis_parameter_list() );    //

                break;

            case moris::Parameter_List_Type::MIG:
                tParameterList = ( moris::prm::create_mig_parameter_list() );

                break;

            case moris::Parameter_List_Type::WRK:
                tParameterList = ( moris::prm::create_wrk_parameter_list() );
                break;

            case moris::Parameter_List_Type::MORISGENERAL:
                switch ( aChild )
                {
                    case 0:
                    {
                        tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                    }
                    break;

                    case 1:
                    {
                        tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                    }
                    break;

                    case 2:
                    {
                        tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                    }
                    break;

                    default:
                        break;
                }
                break;

            default:
                // MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
                break;
        }

        return tParameterList;
    }

}    // namespace moris
#include "main_gui.moc"
