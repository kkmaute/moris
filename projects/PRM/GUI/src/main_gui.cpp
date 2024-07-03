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

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

Moris_Gui::Moris_Gui( QWidget *parent )
        : QWidget( parent )
{

    // mLayout is the main layout of the GUI
    // mSidePanel is the side panel layout where the tree widget and buttons are placed
    mLayout->addLayout( mSidePanel );
    mTreeWidget->setColumnCount( 1 );

    // Project names to be added to mTreeWidget
    mProjectNames = { "OPT", "HMR", "STK", "XTK", "GEN", "FEM", "SOL", "MSI", "VIS", "MIG", "WRK", "MORISGENERAL" };

    // Initialize combo boxes structure
    mComboBox.resize( mProjectNames.size() );

    // Parameter_List_Type ENUM, used to iterate through the projects
    moris::Parameter_List_Type tModule = moris::Parameter_List_Type::OPT;
    QStringList                tStringList;
    QStringList                tSubChildrenList;

    // Initialize data structures for each project
    mTreeWidgetItems.resize( mProjectNames.size() );
    mTreeWidgetChildren.resize( mProjectNames.size() );
    mTreeWidgetSubChildren.resize( mProjectNames.size() );
    mFormLayout.resize( mProjectNames.size() );
    mScrollWidget.resize( mProjectNames.size() );
    mScrollArea.resize( mProjectNames.size() );
    mLineEdit.resize( mProjectNames.size() );
    mCountProps.resize( mProjectNames.size() );


    /* The following for loop loops through the Projects,
    getting the respective outer and inner sub-modules
    and resizing all the lists in that size */

    for ( int iRoot = 0; iRoot < mProjectNames.size(); iRoot++ )
    {

        tStringList = get_outer_sub_parameter_list( tModule );

        // Setting the mTreeWidgetItems with the project names
        mTreeWidgetItems[ iRoot ] = new QTreeWidgetItem;
        mTreeWidgetItems[ iRoot ]->setText( 0, mProjectNames[ iRoot ] );

        // Resizing all the lists for the root index
        mTreeWidgetChildren[ iRoot ].resize( tStringList.size() );
        mTreeWidgetSubChildren[ iRoot ].resize( tStringList.size() );
        mFormLayout[ iRoot ].resize( tStringList.size() );
        mScrollWidget[ iRoot ].resize( tStringList.size() );
        mScrollArea[ iRoot ].resize( tStringList.size() );
        mLineEdit[ iRoot ].resize( tStringList.size() );
        mCountProps[ iRoot ].resize( tStringList.size() );
        mComboBox[ iRoot ].resize( tStringList.size() );

        for ( int iChildren = 0; iChildren < tStringList.size(); iChildren++ )
        {
            /*
            Resizing all the lists for the child index to 1 element for the inner sub-modules.
            Every sub-module only has 1 form but some sub-modules may have more sub-forms (i.e. GEN/Geometries or OPT/Algorithms or SOL/LinearAlgorithms)
            For every sub-module the 0th form will correspond to the main form (i.e. OPT/Algorithms)
                and the rest will correspond to the inner sub-modules (i.e. OPT/Algorithms/gcmma, OPT/Algorithms/lbfgs, etc.
            */
            mTreeWidgetChildren[ iRoot ][ iChildren ] = new QTreeWidgetItem;
            mTreeWidgetChildren[ iRoot ][ iChildren ]->setText( 0, tStringList[ iChildren ] );

            mTreeWidgetSubChildren[ iRoot ][ iChildren ].resize( 1 );
            mFormLayout[ iRoot ][ iChildren ].resize( 1 );
            mScrollWidget[ iRoot ][ iChildren ].resize( 1 );
            mScrollArea[ iRoot ][ iChildren ].resize( 1 );
            mLineEdit[ iRoot ][ iChildren ].resize( 1 );
            mComboBox[ iRoot ][ iChildren ] = new QComboBox();

            /*
            Some sub-modules may not have any inner sub-modules,
            this if statement only appends 1 element
            in each of the lists if that is the case.

            Adding all the items to the combo boxes for the respective sub-module
            */

            QList< int > tCountPropsSubChildList;

            if ( iRoot == static_cast< int >( moris::Parameter_List_Type::OPT ) && iChildren == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) )
            {
                // Adding the OPT/Algorithms to its combo box

                mComboBox[ iRoot ][ iChildren ]->addItems( { "gcmma", "lbfgs", "sql", "sweep" } );
                int tComboBoxCount = mComboBox[ iRoot ][ iChildren ]->count();
                tCountPropsSubChildList.resize( tComboBoxCount + 1 );
            }
            else if ( iRoot == static_cast< int >( moris::Parameter_List_Type::GEN ) && iChildren == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) )
            {
                // Adding the GEN/Geometries to its combo box

                mComboBox[ iRoot ][ iChildren ]->addItems( { "NONE",
                        "CONSTANT",
                        "LINE",
                        "CIRCLE",
                        "SUPERELLIPSE",
                        "PLANE",
                        "SPHERE",
                        "SUPERELLIPSOID",
                        "SCALED_FIELD",
                        "COMBINED_FIELDS",
                        "NODAL",
                        "NODAL_FROM_FILE",
                        "SIGNED_DISTANCE_OBJECT",
                        "SIGNED_DISTANCE_IMAGE",
                        "USER_DEFINED",
                        "SURFACE_MESH",
                        "VOXEL" } );
                int tComboBoxCount = mComboBox[ iRoot ][ iChildren ]->count();
                tCountPropsSubChildList.resize( tComboBoxCount + 1 );
            }
            else if ( iRoot == static_cast< int >( moris::Parameter_List_Type::SOL ) && iChildren == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) )
            {
                // Adding the SOL/LinearAlgorithms to its combo box

                mComboBox[ iRoot ][ iChildren ]->addItems( { "Aztec", "Amesos", "Belos", "PETSC", "EigenSolver", "ML", "Slepc_Solver" } );
                int tComboBoxCount = mComboBox[ iRoot ][ iChildren ]->count();
                tCountPropsSubChildList.resize( tComboBoxCount + 1 );
            }
            else
            {
                // If the sub-module does not have any sub-forms (i.e. gcmma), then only 1 element is added to mCountProps

                tCountPropsSubChildList.resize( 1 );
            }

            // Initializing all the mCountProps to 0
            for ( int &tVal : tCountPropsSubChildList )
            {
                tVal = 0;
            }

            mCountProps[ iRoot ][ iChildren ] = tCountPropsSubChildList;

            setup_scroll_widget( iRoot, iChildren, 0 );

            // Adding elements to the first form layout to show up upon opening the GUI
            if ( iRoot == static_cast< int >( moris::Parameter_List_Type::OPT ) && iChildren == static_cast< int >( moris::OPT_SubModule::OPTIMIZATION_PROBLEMS ) )
            {
                set_form_visible( iRoot, iChildren, 0, true );
                add_elements( 0, 0, 0 );
                mCountProps[ 0 ][ 0 ][ 0 ] = 1;
            }
            else
            {
                set_form_visible( iRoot, iChildren, 0, false );
            }
        }

        // Adding the children to the root items
        mTreeWidgetItems[ iRoot ]->addChildren( mTreeWidgetChildren[ iRoot ] );

        // Incrementing the module to get the next project name
        tModule = static_cast< moris::Parameter_List_Type >( static_cast< int >( tModule ) + 1 );
    }

    // Adding the root items to the tree widget and setting the header label
    mTreeWidget->addTopLevelItems( mTreeWidgetItems );
    mTreeWidget->setHeaderLabel( "Projects" );

    // Adding the tree widget to the side panel
    mSidePanel->addWidget( mTreeWidget );

    // Setting the text for the add and remove buttons
    mAddButton->setText( "Add" );
    mSidePanel->addWidget( mAddButton );

    mRemoveButton->setText( "Remove" );
    mSidePanel->addWidget( mRemoveButton );

    /*
    mOldSelection keeps track of the index of the old selection,
    this is to check if a treeWidgetItem without any form associated (i.e. Project names) is clicked
    */

    int tOldSelection = 0;
    mOldSelection.append( tOldSelection );
    mOldSelection.append( tOldSelection );
    mOldSelection.append( tOldSelection );

    /*
    When the tree widget's current item is changed, the function parameter_selected runs
    When the "Add" button is clicked, the function add_more_props() is ran
    When the "Remove" button is clicked, the function remove_props() is ran
    */

    connect( mTreeWidget, SIGNAL( currentItemChanged( QTreeWidgetItem *, QTreeWidgetItem * ) ), this, SLOT( parameter_selected( QTreeWidgetItem *, QTreeWidgetItem * ) ) );
    connect( mAddButton, SIGNAL( clicked() ), this, SLOT( add_more_props() ) );
    connect( mRemoveButton, SIGNAL( clicked() ), this, SLOT( remove_props() ) );
}

void Moris_Gui::setup_scroll_widget( int aRoot, int aChild, int aSubChild )
{

    /*
    function name: setup_scroll_widget
    parameters:
      int aRoot -> this gives the project name index
      int aChild -> gives the child index
      int aSubChild -> gives the Sub-Child (inner sub-module) index
    returns:
        NONE
    description: This function sets up the scroll functionality for the form layout. Called when a new form is added.
    */

    mScrollWidget[ aRoot ][ aChild ][ aSubChild ] = new QWidget;
    mScrollArea[ aRoot ][ aChild ][ aSubChild ]   = new QScrollArea;
    mFormLayout[ aRoot ][ aChild ][ aSubChild ]   = new QFormLayout;
    mScrollWidget[ aRoot ][ aChild ][ aSubChild ]->setLayout( mFormLayout[ aRoot ][ aChild ][ aSubChild ] );
    mScrollArea[ aRoot ][ aChild ][ aSubChild ]->setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
    mScrollArea[ aRoot ][ aChild ][ aSubChild ]->setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded );
    mScrollArea[ aRoot ][ aChild ][ aSubChild ]->setWidgetResizable( true );
    mScrollArea[ aRoot ][ aChild ][ aSubChild ]->setWidget( mScrollWidget[ aRoot ][ aChild ][ aSubChild ] );
    mLayout->addWidget( mScrollArea[ aRoot ][ aChild ][ aSubChild ] );
}

QStringList Moris_Gui::get_outer_sub_parameter_list( moris::Parameter_List_Type aModule )
{

    /*
    function name: get_outer_sub_parameter_list
    parameters:
      moris::Parameter_List_Type aModule (ENUM) -> this gives the project name
    returns:
        QStringList:
            The list of sub-modules corresponding to the input module is return
    */

    // initialize the names with the standard
    QStringList tNames = { "General" };

    // get the names of the sub-parameter lists for each of the modules
    switch ( aModule )
    {
        case moris::Parameter_List_Type::OPT:
            tNames = { "OptimizationProblems", "Interface", "Algorithms" };
            break;

        case moris::Parameter_List_Type::HMR:
            break;    // standard name

        case moris::Parameter_List_Type::STK:
            break;    // standard name

        case moris::Parameter_List_Type::XTK:
            break;    // standard name

        case moris::Parameter_List_Type::GEN:
            tNames = { "General", "Geometries", "Properties" };
            break;

        case moris::Parameter_List_Type::FEM:
            tNames = {
                "Properties",                 // 0
                "ConstitutiveModels",         // 1
                "StabilizationParameters",    // 2
                "IWG",                        // 3
                "IQI",                        // 4
                "ComputationParameters",      // 5
                "Fields",                     // 6
                //"Materials",                  // 7
                "MaterialModels"    // 8
            };
            break;

        case moris::Parameter_List_Type::SOL:
            tNames = {
                "LinearAlgorithms",        // 0
                "LinearSolvers",           // 1
                "NonLinearAlgorithms",     // 2
                "NonLinearSolvers",        // 3
                "TimeSolverAlgorithms",    // 4
                "TimeSolvers",             // 5
                "SolverWarehouse",         // 6
                "Preconditioners"          // 7
            };
            break;

        case moris::Parameter_List_Type::MSI:
            break;    // standard name

        case moris::Parameter_List_Type::VIS:
            tNames = { "OutputMeshes" };
            break;

        case moris::Parameter_List_Type::MIG:
            break;    // standard name

        case moris::Parameter_List_Type::WRK:
            break;    // standard name

        case moris::Parameter_List_Type::MORISGENERAL:
            // tNames = { "Remeshing", "Refinement", "Mapping" };
            break;

        default:
            // MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
            break;
    }

    // check validity of the input
    // uint tNumSubParamLists = tNames.size();

    // retrieve the name for the specific sub-parameter list requested
    return tNames;
}

void Moris_Gui::set_form_visible( int aRoot, int aChildren, int aSubChildren, bool aCheck )
{
    /*
    function name: set_form_visible
    parameters:
      int aRoot -> this gives the project name index
      int aChild -> gives the child index
      int aSubChild -> gives the Sub-Child (inner sub-module) index
    returns:
        NONE
    */


    for ( int i = 0; i < mFormLayout[ aRoot ][ aChildren ][ aSubChildren ]->rowCount(); i++ )
    {
        mFormLayout[ aRoot ][ aChildren ][ aSubChildren ]->setRowVisible( i, aCheck );
    }
    mScrollWidget[ aRoot ][ aChildren ][ aSubChildren ]->setVisible( aCheck );
    mScrollArea[ aRoot ][ aChildren ][ aSubChildren ]->setVisible( aCheck );
}

QList< QStringList > Moris_Gui::get_parameter_list( moris::Parameter_List_Type aModule, int aChild, int aSubChild )
{
    /*
    function name: get_parameter_list
    parameters:
      moris::Parameter_List_Type aModule (ENUM) -> this gives the project name
      int aChild -> gives the child index
      int aSubChild -> gives the Sub-Child (inner sub-module) index
    returns:
        QList <QStringList>
            the create_function returns a ParameterList object that is a type of map
            The 0th index of the QList gives the "keys" of the map
            The 1st index of the QList gives the default "values" of the map
    */


    QList< QStringList > tParameterList;

    switch ( aModule )
    {
        case moris::Parameter_List_Type::OPT:
            switch ( aChild )
            {
                case 0:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_opt_problem_parameter_list() );
                    break;

                case 1:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_opt_interface_parameter_list() );

                    // Commented out the Interface manager for now

                    // switch ( aSubChild )
                    // {
                    //     case 0:
                    //         tParameterList = convert_parameters_to_QStringList( moris::prm::create_opt_interface_parameter_list() );

                    //         break;

                    //     case 1:
                    //         tParameterList = convert_parameters_to_QStringList( moris::prm::create_opt_interface_manager_parameter_list() );

                    //         break;

                    //     default:
                    //         break;
                    // }

                    break;

                case 2:
                    switch ( aSubChild )
                    {
                        case 1:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_gcmma_parameter_list() );

                            break;

                        case 2:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_lbfgs_parameter_list() );
                            break;

                        case 3:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_sqp_parameter_list() );

                            break;

                        case 4:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_sweep_parameter_list() );

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
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_hmr_parameter_list() );

            break;

        case moris::Parameter_List_Type::STK:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_stk_parameter_list() );
            break;

        case moris::Parameter_List_Type::XTK:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_xtk_parameter_list() );
            break;

        case moris::Parameter_List_Type::GEN:
            switch ( aChild )
            {
                case 0:
                {
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_gen_parameter_list() );
                    break;
                }

                case 1:
                {
                    int aRoot = static_cast< int >( aModule );
                    if ( mComboBox[ aRoot ][ aChild ]->currentIndex() <= 14 )
                    {
                        tParameterList = convert_parameters_to_QStringList( moris::prm::create_level_set_geometry_parameter_list( static_cast< moris::gen::Field_Type >( mComboBox[ aRoot ][ aChild ]->currentIndex() ) ) );
                    }
                    else if ( mComboBox[ aRoot ][ aChild ]->currentIndex() == 15 )
                    {
                        tParameterList = convert_parameters_to_QStringList( moris::prm::create_surface_mesh_geometry_parameter_list() );
                    }
                    else
                    {
                        tParameterList = convert_parameters_to_QStringList( moris::prm::create_voxel_geometry_parameter_list() );
                    }
                    break;
                }
                case 2:
                {
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_gen_property_parameter_list( moris::gen::Field_Type::CONSTANT ) );
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
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_property_parameter_list() );
                    break;

                case 1:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_constitutive_model_parameter_list() );
                    break;

                case 2:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_stabilization_parameter_parameter_list() );
                    break;

                case 3:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_IWG_parameter_list() );
                    break;

                case 4:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_IQI_parameter_list() );
                    break;

                case 5:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_computation_parameter_list() );
                    break;

                case 6:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_fem_field_parameter_list() );
                    break;

                case 7:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_material_model_parameter_list() );
                    break;

                default:
                    break;
            }


            break;

        case moris::Parameter_List_Type::SOL:
            tParameterList.resize( 8 );

            switch ( aChild )
            {
                case 0:

                    switch ( aSubChild )
                    {
                        case 0:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_algorithm_parameter_list_aztec() );
                            break;

                        case 1:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_algorithm_parameter_list_amesos() );
                            break;

                        case 2:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_algorithm_parameter_list_belos() );
                            break;

                        case 3:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_algorithm_parameter_list_petsc() );
                            break;

                        case 4:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_eigen_algorithm_parameter_list() );
                            break;

                        case 5:
                            // Need to add ML here
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_algorithm_parameter_list_belos() );
                            break;

                        case 6:
                            tParameterList = convert_parameters_to_QStringList( moris::prm::create_slepc_algorithm_parameter_list() );
                            break;

                        default:
                            break;
                    }
                    break;

                    break;

                case 1:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_linear_solver_parameter_list() );
                    break;

                case 2:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_nonlinear_algorithm_parameter_list() );
                    break;

                case 3:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_nonlinear_solver_parameter_list() );
                    break;

                case 4:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_time_solver_algorithm_parameter_list() );
                    break;

                case 5:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_time_solver_parameter_list() );
                    break;

                case 6:
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_solver_warehouse_parameterlist() );
                    break;

                case 7:
                    // Need to add Preconditioners
                    tParameterList = convert_parameters_to_QStringList( moris::prm::create_material_model_parameter_list() );
                    break;

                default:
                    break;
            }

            break;

        case moris::Parameter_List_Type::MSI:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_msi_parameter_list() );
            break;

        case moris::Parameter_List_Type::VIS:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_vis_parameter_list() );    //

            break;

        case moris::Parameter_List_Type::MIG:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_mig_parameter_list() );

            break;

        case moris::Parameter_List_Type::WRK:
            tParameterList = convert_parameters_to_QStringList( moris::prm::create_wrk_parameter_list() );
            break;

        case moris::Parameter_List_Type::MORISGENERAL:
            break;

        default:
            // MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
            break;
    }

    return tParameterList;
}

void Moris_Gui::add_elements( int aRoot, int aChild, int aSubChild )
{

    /*
    function name: add_elements
    parameters:
      int aRoot -> this gives the project name index
      int aChild -> gives the child index
      int aSubChild -> gives the Sub-Child (inner sub-module) index
    returns:
        NONE
    description: This function actually adds elements like when given the index where to put
    */

    /*
    The following if statement checks if the index is for OPT/Algorithm or GEN/Geometries or SOL/LinearAlgorithms.
     If it is, then it adds a combo box with the respective options,
     if not then it adds the respective QLineEdit fields by calling get_parameter_list.
    */
    if ( ( ( aRoot == static_cast< int >( moris::Parameter_List_Type::OPT ) && aChild == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) && aSubChild == 0 )
                 || ( aRoot == static_cast< int >( moris::Parameter_List_Type::GEN ) && aChild == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) && aSubChild == 0 )
                 || ( aRoot == static_cast< int >( moris::Parameter_List_Type::SOL ) && aChild == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) && aSubChild == 0 ) )
            && mCountProps[ aRoot ][ aChild ][ aSubChild ] == 0 )
    {
        // Adding a row with the combo box to the form layout for OPT/Algorithm, GEN/Geometries, SOL/LinearAlgorithms
        mFormLayout[ aRoot ][ aChild ][ aSubChild ]->addRow( "Please select:", mComboBox[ aRoot ][ aChild ] );

        // mCountProps of the respective index is incremented so that this is not called again
        mCountProps[ aRoot ][ aChild ][ aSubChild ]++;
    }
    else
    {
        // Get the parameter list for the respective module, child and sub-child

        moris::Parameter_List_Type aModule     = static_cast< moris::Parameter_List_Type >( aRoot );
        QList< QStringList >       tStringList = get_parameter_list( aModule, aChild, aSubChild );

        // Check the number of elements already present in mLineEdit
        int tCounter = mLineEdit[ aRoot ][ aChild ][ aSubChild ].size();

        // If the index is for OPT/Algorithm or GEN/Geometries or SOL/LinearAlgorithms, then add a new sub-child form associated elements for the layout
        // and then add the QLineEdit fields to the sub-form layout

        // If not then only the QLineEdit fields are added to the respective form layout
        if ( ( aRoot == static_cast< int >( moris::Parameter_List_Type::OPT ) && aChild == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) )
                || ( aRoot == static_cast< int >( moris::Parameter_List_Type::GEN ) && aChild == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) )
                || ( aRoot == static_cast< int >( moris::Parameter_List_Type::SOL ) && aChild == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) ) )
        {
            /*
            In ths case we are in the OPT/Algorithm, GEN/Geometries, SOL/LinearAlgorithms sub-modules and adding a new form base on the combo box selection
            Thus new forms, layouts, tree widgets, scroll widgets, and line edits are added.
            */
            mFormLayout[ aRoot ][ aChild ].append( new QFormLayout );

            mTreeWidgetSubChildren[ aRoot ][ aChild ].append( new QTreeWidgetItem );

            mScrollWidget[ aRoot ][ aChild ].append( new QWidget );
            mScrollArea[ aRoot ][ aChild ].append( new QScrollArea );

            QList< QLineEdit * > tLineEditChild;
            mLineEdit[ aRoot ][ aChild ].append( tLineEditChild );

            // aSubChild is set to the last index of the sub-child list after we add the new forms
            aSubChild = mFormLayout[ aRoot ][ aChild ].size() - 1;

            // Setting up the scroll widget for the new form
            setup_scroll_widget( aRoot, aChild, aSubChild );

            /*
            The following few lines are setting the name of the new form in the tree widget based on the combo box selection
            and the number of those forms already present.
            Then adding that form to the tree widget and adding the form layout to the main layout
            */

            QString tParameterName = mComboBox[ aRoot ][ aChild ]->currentText() + QString::number( mCountProps[ aRoot ][ aChild ][ mComboBox[ aRoot ][ aChild ]->currentIndex() + 1 ] );
            mTreeWidgetSubChildren[ aRoot ][ aChild ][ aSubChild ]->setText( 0, tParameterName );
            mCountProps[ aRoot ][ aChild ][ mComboBox[ aRoot ][ aChild ]->currentIndex() + 1 ]++;

            mTreeWidgetChildren[ aRoot ][ aChild ]->addChild( mTreeWidgetSubChildren[ aRoot ][ aChild ][ aSubChild ] );

            mLayout->addWidget( mScrollArea[ aRoot ][ aChild ][ aSubChild ] );

            set_form_visible( aRoot, aChild, aSubChild, false );

            tCounter = mLineEdit[ aRoot ][ aChild ][ aSubChild ].size();
        }
        else
        {
            // If it is not in the "special" sub-modules then it is one of the general cases and the mCountProps of that element in added by 1
            mCountProps[ aRoot ][ aChild ][ aSubChild ]++;
        }

        for ( int iElements = 0; iElements < tStringList[ 0 ].size(); iElements++ )
        {
            // Create a temp lineEdit field, set placeholder to that and then append to the list
            mLineEdit[ aRoot ][ aChild ][ aSubChild ].append( new QLineEdit() );

            // Only setting the place holder text if it is not an empty string
            if ( tStringList[ 1 ][ iElements ].toStdString() != "\"\"" )
            {
                mLineEdit[ aRoot ][ aChild ][ aSubChild ][ iElements + tCounter ]->setPlaceholderText( tStringList[ 1 ][ iElements ] );
            }

            // Adding the line edit field to the form layout,
            // if it is one of the "special" cases the it is added to the last form layout
            // if it is one of the general cases then it is added to the current form layout selected
            if ( ( aRoot == 0 && aChild == 2 ) || ( aRoot == 4 && aChild == 1 ) || ( aRoot == 6 && aChild == 0 ) )
            {
                mFormLayout[ aRoot ][ aChild ][ mFormLayout[ aRoot ][ aChild ].size() - 1 ]->addRow( tStringList[ 0 ][ iElements ], mLineEdit[ aRoot ][ aChild ][ aSubChild ][ iElements + tCounter ] );
            }
            else
            {
                mFormLayout[ aRoot ][ aChild ][ aSubChild ]->addRow( tStringList[ 0 ][ iElements ], mLineEdit[ aRoot ][ aChild ][ aSubChild ][ iElements + tCounter ] );
            }
        }
    }

    // mCountProps[ aRoot ][ aChild ][ aSubChild ]++;
}

void Moris_Gui::parameter_selected( QTreeWidgetItem *aNewItem, QTreeWidgetItem *aOldItem )
{
    /*
    function name: parameter_selected
    parameters:
      aNewItem (QTreeWidgetItem) -> the new selection by the user
      aOldItem (QTreeWidgetItem)-> the old selection by the user
    returns:
        NONE
    description: This function is the SLOT for SIGNAL from currentItemChanged( QTreeWidgetItem *, QTreeWidgetItem * ) in mTreeWidget
    It checks the new and old indices by calling get_tree_index, checks if there is a corresponding layout for the previous/new item and calls the mOldSelection as needed
    */

    QList< int > tNewIndex = get_tree_index( aNewItem );
    QList< int > tOldIndex;

    if ( aOldItem )
    {
        tOldIndex = get_tree_index( aOldItem );
    }
    else
    {
        return;
    }

    if ( tOldIndex.isEmpty() && tNewIndex.isEmpty() )
    {
        return;
    }


    // If the given new index is not empty and the number of properties is not 0 (i.e. this is the first time this form is being selected)
    // and the form is not one of the sub-forms of the special cases (i.e. OPT/Algorithms, GEN/Geometries, SOL/LinearAlgorithms)
    // then add the elements to the form.
    if ( !tNewIndex.isEmpty()
            && mCountProps[ tNewIndex[ 0 ] ][ tNewIndex[ 1 ] ][ tNewIndex[ 2 ] ] == 0
            && !( ( tNewIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::OPT ) && tNewIndex[ 1 ] == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) && tNewIndex[ 2 ] != 0 )
                    || ( tNewIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::GEN ) && tNewIndex[ 1 ] == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) && tNewIndex[ 2 ] != 0 )
                    || ( tNewIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::SOL ) && tNewIndex[ 1 ] == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) && tNewIndex[ 2 ] != 0 ) ) )
    {
        add_elements( tNewIndex[ 0 ], tNewIndex[ 1 ], tNewIndex[ 2 ] );
    }

    int tNewRoot = 0, tNewChild = 0, tNewSubChild = 0;
    int tOldRoot = 0, tOldChild = 0, tOldSubChild = 0;

    if ( !tOldIndex.isEmpty() )
    {
        tOldRoot     = tOldIndex[ 0 ];
        tOldChild    = tOldIndex[ 1 ];
        tOldSubChild = tOldIndex[ 2 ];
    }
    else
    {
        set_form_visible( mOldSelection[ 0 ], mOldSelection[ 1 ], mOldSelection[ 2 ], false );
    }

    if ( !tNewIndex.isEmpty() )
    {
        tNewRoot     = tNewIndex[ 0 ];
        tNewChild    = tNewIndex[ 1 ];
        tNewSubChild = tNewIndex[ 2 ];

        // Saving the selection for the next time if the user next selects a treeWidgetItem without any form associated (i.e. Project names)
        mOldSelection = tNewIndex;
    }
    else
    {
        // qDebug() << "No matching index found.";
        return;
    }

    set_form_visible( tOldRoot, tOldChild, tOldSubChild, false );
    set_form_visible( tNewRoot, tNewChild, tNewSubChild, true );
}

QList< int > Moris_Gui::get_tree_index( QTreeWidgetItem *tItem )
{

    /*
    function name: get_tree_index
    parameters:
      tItem (QTreeWidgetItem) -> the QTreeWidgetItem selected from mTreeWidget
    returns:
        QList <int> tIndex;
    description: This function returns the 3D index (format: root, sub-child, inner sub-child) of the item passed into the function
        tIndex[0] -> gives the root
        tIndex[1] -> gives the sub-child
        tIndex[2] -> gives the inner sub-child
    */

    QList< int > tIndex;

    // Checking across the root and children items. If the item is still not found then checking the sub-children
    for ( int iRoot = 0; iRoot < mTreeWidgetChildren.size(); iRoot++ )
    {
        for ( int iChildren = 0; iChildren < mTreeWidgetChildren[ iRoot ].size(); iChildren++ )
        {
            // Check if the current item is the direct child
            if ( tItem == mTreeWidgetChildren[ iRoot ][ iChildren ] )
            {
                tIndex.append( iRoot );
                tIndex.append( iChildren );
                tIndex.append( 0 );    // No subChild in this case
                return tIndex;
            }

            // Check if there are subChildren and look for the item among them
            for ( int iSubChildren = 0; iSubChildren < mTreeWidgetSubChildren[ iRoot ][ iChildren ].size(); iSubChildren++ )
            {
                if ( tItem == mTreeWidgetSubChildren[ iRoot ][ iChildren ][ iSubChildren ] )
                {
                    tIndex.append( iRoot );
                    tIndex.append( iChildren );
                    tIndex.append( iSubChildren );
                    return tIndex;
                }
            }
        }
    }


    qDebug() << "No matching index found.";

    return tIndex;
}

QList< QStringList > Moris_Gui::convert_parameters_to_QStringList( moris::Parameter_List tList )
{

    /*
    function name: convert_parameters_to_QStringList
    parameters:
      tList (moris::Parameter_List) -> takes in the Parameter_List map object from MORIS
    returns:
        QList <QStringList> tQtList;
    description: takes in the Parameter_List object from MORIS and converts it into a QList <QStringList> where
        tQtList[0] corresponds to the "key" values of the parameter map
        tQtList[1] corresponds to the "value" of the parameter map
    */

    QList< QStringList > tQtList;
    QStringList          tStringList;
    tQtList.append( tStringList );
    tQtList.append( tStringList );

    for ( auto it = tList.begin(); it != tList.end(); ++it )
    {
        tQtList[ 0 ].append( QString::fromStdString( it->first ) );
        tQtList[ 1 ].append( QString::fromStdString( it->second.get_string() ) );
    }

    return tQtList;
}

void Moris_Gui::add_more_props()
{

    /*
    function name: add_more_props
    parameters:
      NONE
    returns:
        NONE
    description: This is the SLOT for the SIGNAL clicked() from the OBJECT mAddButton
        This function runs whenever the mAddButton is clicked,
        the function checks the index and then depending on where we are in the index it runs add_elements
    */

    QList< int > tIndex = get_tree_index( mTreeWidget->currentItem() );

    if ( ( tIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::OPT ) && tIndex[ 1 ] == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) )
            || ( tIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::GEN ) && tIndex[ 1 ] == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) )
            || ( tIndex[ 0 ] == static_cast< int >( moris::Parameter_List_Type::SOL ) && tIndex[ 1 ] == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) ) )
    {
        // If we are in the special cases then we need to add a new form based on the combo box selection, so the tIndex[2] is set to the combo box index
        tIndex[ 2 ] = mComboBox[ tIndex[ 0 ] ][ tIndex[ 1 ] ]->currentIndex() + 1;

        add_elements( tIndex[ 0 ], tIndex[ 1 ], tIndex[ 2 ] );
    }
    else
    {
        // General case for adding elements given the index
        add_elements( tIndex[ 0 ], tIndex[ 1 ], tIndex[ 2 ] );
    }

    // mCountProps[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ tIndex[ 2 ] ]++;
}

void Moris_Gui::remove_props()
{

    QList< int > tIndex = get_tree_index( mTreeWidget->currentItem() );

    if ( tIndex.isEmpty() || tIndex.size() < 3 )
    {
        QMessageBox::warning( this, "Error", "Invalid selection." );
        return;
    }

    int tRoot     = tIndex[ 0 ];
    int tChild    = tIndex[ 1 ];
    int tSubChild = tIndex[ 2 ];

    qDebug() << "subChild index" << tSubChild;

    // Extract the name of the current item and remove the last character
    QString tCurrentItemName = mTreeWidget->currentItem()->text( 0 );
    if ( !tCurrentItemName.isEmpty() )
    {
        tCurrentItemName.chop( 1 );    // Remove the last character
    }

    int tCorrectSubChild = tSubChild;

    // Special cases for OPT/Algorithm, GEN/Geometries, SOL/LinearAlgorithms
    if ( ( tRoot == static_cast< int >( moris::Parameter_List_Type::OPT ) && tChild == static_cast< int >( moris::OPT_SubModule::ALGORITHMS ) )
            || ( tRoot == static_cast< int >( moris::Parameter_List_Type::GEN ) && tChild == static_cast< int >( moris::GEN_SubModule::GEOMETRIES ) )
            || ( tRoot == static_cast< int >( moris::Parameter_List_Type::SOL ) && tChild == static_cast< int >( moris::SOL_SubModule::LINEAR_ALGORITHMS ) ) )
    {
        if ( tSubChild == 0 )
        {
            // If the user is trying to remove the first form (i.e. Algorithm form where the combo box is present) then this prevents it
            QMessageBox::warning( this, "Error", "Cannot remove this form" );
            return;
        }

        // Finding the index of the form to be removed in the combo box to keep track of number of elements
        // tSubChild is the form to be deleted and tCorrectSubChild is the index of the item in the combo box to modify mCountProps

        for ( int i = 0; i < mComboBox[ tRoot ][ tChild ]->count(); ++i )
        {
            if ( mComboBox[ tRoot ][ tChild ]->itemText( i ) == tCurrentItemName )
            {
                tCorrectSubChild = i + 1;
                break;
            }
        }

        // If the user is trying to remove the last form then select the previous form
        if ( tSubChild == mTreeWidgetSubChildren[ tRoot ][ tChild ].size() - 1 )
        {
            parameter_selected( mTreeWidgetSubChildren[ tRoot ][ tChild ][ tSubChild - 1 ], mTreeWidget->currentItem() );
        }

        // Remove the form elements, layout, tree widget, scroll widget, and line edits
        for ( QLineEdit *lineEdit : mLineEdit[ tRoot ][ tChild ][ tSubChild ] )
        {
            delete lineEdit;
        }
        mLineEdit[ tRoot ][ tChild ].removeAt( tSubChild );

        delete mFormLayout[ tRoot ][ tChild ].takeAt( tSubChild );
        delete mScrollWidget[ tRoot ][ tChild ].takeAt( tSubChild );
        delete mTreeWidgetSubChildren[ tRoot ][ tChild ].takeAt( tSubChild );
        delete mScrollArea[ tRoot ][ tChild ].takeAt( tSubChild );

        // Decrement the count of properties for the item in the combo box
        mCountProps[ tRoot ][ tChild ][ tCorrectSubChild ]--;

        // Make sure the indices are still valid after removal
        if ( tSubChild > 0 )
        {
            QList< int > tNewIndex = get_tree_index( mTreeWidget->currentItem() );
            set_form_visible( tNewIndex[ 0 ], tNewIndex[ 1 ], tNewIndex[ 2 ], true );
        }
    }
    else
    {
        // General case for removing elements
        if ( mCountProps[ tRoot ][ tChild ][ tSubChild ] <= 1 )
        {
            QMessageBox::warning( this, "Warning", "Cannot remove more properties." );
            return;
        }

        int tIterator = mLineEdit[ tRoot ][ tChild ][ tSubChild ].size() / mCountProps[ tRoot ][ tChild ][ tSubChild ];
        for ( int iRemove = 0; iRemove < tIterator; iRemove++ )
        {
            mFormLayout[ tRoot ][ tChild ][ tSubChild ]->removeRow( mFormLayout[ tRoot ][ tChild ][ tSubChild ]->rowCount() - 1 );
            mLineEdit[ tRoot ][ tChild ][ tSubChild ].removeLast();
        }

        mCountProps[ tRoot ][ tChild ][ tSubChild ]--;
    }
}


#include "main_gui.moc"
