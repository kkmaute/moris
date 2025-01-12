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


// Comm_Manager gMorisComm;
// Logger       gLogger;
namespace moris
{
    // uint       tag       = 0;
    // Moris_Gui *gMorisGui = nullptr;
    Moris_Gui::Moris_Gui( QWidget *parent )
            : QWidget( parent )
    {
        this->setWindowTitle( "MORIS GUI" );
        this->resize( 800, 600 );
        // gMorisGui         = this;
        QString tFilePath = get_moris_file_path();

        if ( tFilePath.isEmpty() )
        {
            mLibrary.create_new_module_parameterlist();
        }
        else
        {
            // load the parameter list from the xml file
            mLibrary.load_parameter_list( tFilePath.toStdString(), File_Type::XML_FILE );

            // Update mParameterLists (sitting in cl_Library_IO) with the parameter lists from the xml file
            mLibrary.load_parameters_from_xml();
        }

        // Initialize the GUI
        initialize_gui();
    }

    Moris_Gui::Moris_Gui( QWidget *parent, std::string aFileName )
            : QWidget( parent )
    {
        // gMorisGui         = this;
        QString tFilePath = QString::fromStdString( aFileName );

        if ( tFilePath.isEmpty() )
        {
            mLibrary.create_new_module_parameterlist();
        }
        else
        {
            // load the parameter list from the xml file
            mLibrary.load_parameter_list( tFilePath.toStdString(), File_Type::XML_FILE );

            // Update mParameterLists (sitting in cl_Library_IO) with the parameter lists from the xml file
            mLibrary.load_parameters_from_xml();
        }

        // Initialize the GUI
        initialize_gui();
    }

    void Moris_Gui::initialize_gui()
    {

        //  mLayout is the main layout of the GUI
        //  mSidePanel is the side panel layout where the tree widget and buttons are placed

        mLayout->addLayout( mSidePanel );
        mTreeWidget->setColumnCount( 1 );

        // Read from mParamterLists instead of moris functions, put these functions in the read() function

        for ( uint iProject = 0; iProject < (uint)( Module_Type::END_ENUM ); iProject++ )
        {
            mProjectNames.append( QString::fromStdString( convert_parameter_list_enum_to_string( (Module_Type)( iProject ) ) ) );
        }

        Module_Type tModule = Module_Type::OPT;
        QStringList tStringList;
        QStringList tSubChildrenList;

        // Initialize data structures for each project
        mTreeWidgetItems.resize( mProjectNames.size() );
        mTreeWidgetChildren.resize( mProjectNames.size() );
        mTreeWidgetSubChildren.resize( mProjectNames.size() );

        // QTreeWidgetItems
        mQTreeWidgetItems.resize( mProjectNames.size() );
        mQTreeWidgetChildren.resize( mProjectNames.size() );
        mQTreeWidgetSubChildren.resize( mProjectNames.size() );


        for ( uint iRoot = 0; iRoot < mLibrary.get_parameter_lists().size(); iRoot++ )
        {
            // Get the number of sub-parameter lists in the module
            tStringList.clear();

            // Read from mParamterLists instead of moris functions, put these functions in the read() function
            Vector< std::string > tSubmoduleNames = get_submodule_names( tModule );
            for ( uint iList = 0; iList < get_number_of_sub_parameter_lists_in_module( tModule ); iList++ )
            {
                tStringList.append( QString::fromStdString( tSubmoduleNames( iList ) ) );
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

                if ( iRoot == (uint)( Module_Type::OPT ) && iChildren == (uint)( OPT_Submodule::OPTIMIZATION_PROBLEMS ) )
                {
                    // Setting the OPT/Optimization Problems form visible upon construction of the GUI
                    // mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mLibrary.get_parameter_lists()( iRoot )(iChildren)( 0 ) );
                    // mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
                    // mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( true );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    mOldItem = mQTreeWidgetChildren[ iRoot ][ iChildren ];
                }
                else if ( iRoot == (uint)( Module_Type::OPT ) && iChildren == (uint)( OPT_Submodule::ALGORITHMS ) )
                {
                    // Adding the OPT/Algorithms to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( { "gcmma", "lbfgs", "sql", "sweep" } );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );

                    // If there are more than one algorithms from the xml file, add the respective forms to the GUI
                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                }
                else if ( iRoot == (uint)( Module_Type::GEN ) && iChildren == (uint)( GEN_Submodule::GEOMETRIES ) )
                {
                    // Adding the GEN/Geometries to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    QStringList tGENGeometriesList;
                    for ( uint iGeometries = 0; iGeometries < (uint)gen::Geometry_Type_String::values.size(); iGeometries++ )
                    {
                        // For Level Set geometries, iterates through the Field_Types
                        if ( (uint)iGeometries == (uint)gen::Geometry_Type::LEVEL_SET )
                        {
                            for ( uint iLevelSet = 0; iLevelSet < (uint)gen::Field_Type_String::values.size(); iLevelSet++ )
                            {
                                tGENGeometriesList.append( QString::fromStdString( gen::Field_Type_String::values( iLevelSet ) ) );
                            }
                        }
                        else
                        {
                            // For other geometries, iterates through the Geometry_Types
                            tGENGeometriesList.append( QString::fromStdString( gen::Geometry_Type_String::values( iGeometries ) ) );
                        }
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tGENGeometriesList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );

                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                }
                else if ( iRoot == (uint)( Module_Type::GEN ) && iChildren == (uint)( GEN_Submodule::PROPERTIES ) )
                {
                    // Adding the GEN/Field to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    QStringList tGENFieldList;
                    for ( uint iField = 0; iField < (uint)gen::Field_Type_String::values.size(); iField++ )
                    {
                        tGENFieldList.append( QString::fromStdString( gen::Field_Type_String::values( iField ) ) );
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tGENFieldList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );

                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                }
                else if ( iRoot == (uint)( Module_Type::SOL ) && iChildren == (uint)( SOL_Submodule::LINEAR_ALGORITHMS ) )
                {
                    // Adding the SOL/Linear Algorithms to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    QStringList tSolverList;
                    for ( uint iSolver = 0; iSolver < (uint)sol::SolverType_String::values.size(); iSolver++ )
                    {
                        tSolverList.append( QString::fromStdString( sol::SolverType_String::values( iSolver ) ) );
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tSolverList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );

                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                }
                else if ( iRoot == (uint)( Module_Type::SOL ) && iChildren == (uint)( SOL_Submodule::PRECONDITIONERS ) )
                {
                    // Adding the SOL/Preconditioners to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    QStringList tPreconditionerList;
                    for ( uint iPreconditioner = 0; iPreconditioner < (uint)sol::PreconditionerType_String::values.size(); iPreconditioner++ )
                    {
                        tPreconditionerList.append( QString::fromStdString( sol::PreconditionerType_String::values( iPreconditioner ) ) );
                    }
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setComboBoxItems( tPreconditionerList );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( false );

                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                }
                else
                {
                    if ( iChildren < mLibrary.get_parameter_lists()( iRoot ).size() )
                    {
                        if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() >= 1 )
                        {
                            for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                            {
                                add_project( iRoot, iChildren, i );
                            }
                        }
                        else
                        {
                            mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                        }
                    }
                    else
                    {
                        mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                        // if ( iRoot == (uint)( Module_Type::FEM ) )
                        // {
                        //     mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                        // }
                        // else if ( iRoot == (uint)( Module_Type::OPT ) && iChildren == (uint)( OPT_Submodule::INTERFACE ) )
                        // {
                        //     mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                        // }
                        // else
                        // {
                        //     // Adding the parameter_list elements to all the forms and setting the visibility to false
                        //     mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mLibrary.get_parameter_lists()( iRoot )(iChildren)( 0 ) );
                        //     mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
                        // }
                    }
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

            tModule = (Module_Type)( (uint)( tModule ) + 1 );

            // MORIS_ERROR( false, "error here" );
        }

        // Setting the header label for the tree widget
        mTreeWidget->setHeaderLabel( "Projects" );

        // Adding the tree widget to the side panel
        mSidePanel->addWidget( mTreeWidget );


        // // Create an instance of Library_IO_Standard
        // Library_IO_Standard library;
        // // Assign the mParameterLists from main_gui to the library instance
        // library.mParameterLists = main_gui::mParameterLists;

        // Setting the text for the add and remove buttons
        mAddButton->setText( "Add" );
        mSidePanel->addWidget( mAddButton );

        mRemoveButton->setText( "Remove" );
        mSidePanel->addWidget( mRemoveButton );

        // Setting the text for the Save to XML button
        m_save_to_xml_Button->setText( "Save to XML" );
        mSidePanel->addWidget( m_save_to_xml_Button );

        // Connecting the signals and slots
        connect( mTreeWidget, SIGNAL( itemSelectionChanged() ), this, SLOT( parameter_selected() ) );
        connect( mAddButton, SIGNAL( clicked() ), this, SLOT( add_more_props() ) );
        connect( mRemoveButton, SIGNAL( clicked() ), this, SLOT( remove_props() ) );
        connect( m_save_to_xml_Button, SIGNAL( clicked() ), this, SLOT( write_to_xml() ) );
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
                uint          tRoot         = tCurrentIndex[ 0 ];
                uint          tChild        = tCurrentIndex[ 1 ];

                // Append a new QTreeWidgetItem and Moris_Tree_Widget_Item to the sub-children list of the current widget
                // And set up the necessary elements

                mQTreeWidgetSubChildren[ tRoot ][ tChild ].append( new QTreeWidgetItem( mQTreeWidgetChildren[ tRoot ][ tChild ] ) );
                mTreeWidgetSubChildren[ tRoot ][ tChild ].append( new Moris_Tree_Widget_Item() );

                mQTreeWidgetChildren[ tRoot ][ tChild ]->addChild( mQTreeWidgetSubChildren[ tRoot ][ tChild ].last() );

                // Reorganizing the count of the sub-forms from the item selection in the combo box
                // All these depend on the item selection in the combo box
                mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->addSubFormCountProps();
                if ( mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->getComboBox()->count() == 0 )
                {
                    // For ordinary sub-forms
                    std::string tSubParameterListName = get_inner_sub_parameter_list_name( (Module_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ] );
                    mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setText( 0, QString::fromStdString( tSubParameterListName ) + " " + QString::number( tCurrentWidget->getSubFormCountProps() ) );
                    mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setSubFormType( -1 );
                    mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).add_parameter_list( create_parameter_list( (Module_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], 0 ) );
                }
                else
                {
                    // For special sub-forms (e.g. OPT/Algorithms, GEN/Geometries, GEN/Field, SOL/Linear Algorithms, SOL/Preconditioners)
                    mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setText( 0, tCurrentWidget->getComboBox()->currentText() + " " + QString::number( tCurrentWidget->getSubFormCountProps() ) );
                    mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setSubFormType( tCurrentWidget->getComboBox()->currentIndex() );
                    // mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).add_parameter_list( create_parameter_list( (Module_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], tCurrentWidget->getComboBox()->currentIndex() ) );
                    mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).add_parameter_list();
                }

                // Setting up the scroll area, index, saving type for the sub-forms
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setupScrollArea();
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setIndex( { tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], (uint)mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].size() - 1 } );

                // Append to mLibrary.get_parameter_lists() and pass it into add_elements
                Submodule_Parameter_Lists &tSubmoduleParameterList = mLibrary.get_parameter_lists()( tRoot )( tChild );
                mTreeWidgetSubChildren[ tRoot ][ tChild ].last()->add_elements( tSubmoduleParameterList( tSubmoduleParameterList.size() - 1 ) );
                mTreeWidgetSubChildren[ tRoot ][ tChild ].last()->set_form_visible( false );

                // Connecting the QTreeWidgetItem to the Moris_Tree_Widget_Item
                mTreeWidget->setItemWidget( mQTreeWidgetSubChildren[ tRoot ][ tChild ].last(), 0, mTreeWidgetSubChildren[ tRoot ][ tChild ].last() );

                // Adding the scroll area to the main layout
                mLayout->addWidget( mTreeWidgetSubChildren[ tRoot ][ tChild ].last()->getScrollArea() );

                if ( tRoot == (uint)Module_Type::FEM && tChild == (uint)FEM_Submodule::PROPERTIES && mQTreeWidgetSubChildren[ tRoot ][ tChild ].size() > 0 )
                {
                    mPropertyNameList.resize( mQTreeWidgetSubChildren[ tRoot ][ tChild ].size() );
                    for ( uint i = 0; i < mQTreeWidgetSubChildren[ tRoot ][ tChild ].size(); i++ )
                    {
                        mPropertyNameList[ i ] = mQTreeWidgetSubChildren[ tRoot ][ tChild ][ i ]->text( 0 );
                    }
                }

                if ( mTreeWidgetItems.size() >= (uint)Module_Type::FEM )
                {
                    if ( mTreeWidgetChildren[ (uint)Module_Type::FEM ].size() >= (uint)FEM_Submodule::CONSTITUTIVE_MODELS )
                    {
                        if ( mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size() > 0 )
                        {
                            for ( uint i = 0; i < mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size(); i++ )
                            {
                                if ( !mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ][ i ]->isPropertyListSet() )
                                {
                                    mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ][ i ]->setPropertyNameList( mPropertyNameList );
                                }
                                // pass the updated property name list to the Moris_Tree_Widget_Item by reference
                            }
                        }
                    }
                }

                // Get mWidget from the Moris_Tree_Widget_Item by reference
                // Loop over all the widgets in mTreeWidgetSubChildren[ tRoot ][ tChild ] and assign tWidget to every one
                QList< QList< QWidget * > > tWidget;
                tWidget.resize( mTreeWidgetSubChildren[ tRoot ][ tChild ].size() );
                for ( uint i = 0; i < mTreeWidgetSubChildren[ tRoot ][ tChild ].size(); i++ )
                {
                    tWidget[ i ] = mTreeWidgetSubChildren[ tRoot ][ tChild ][ i ]->mWidget;

                    for ( uint it = 0; it < tWidget[ i ].size(); it++ )
                    {
                        if ( endswith( tWidget[ i ][ it ]->objectName().toStdString(), "name" ) )
                        {
                            auto &tLineEdit      = dynamic_cast< Moris_Line_Edit      &>( *tWidget[ i ][ it ] );
                            auto &treeWidgetItem = mTreeWidgetSubChildren[ tRoot ][ tChild ][ i ];

                            if ( tRoot == (uint)Module_Type::FEM && tChild == (uint)FEM_Submodule::PROPERTIES && mQTreeWidgetSubChildren[ tRoot ][ tChild ].size() > 0 )
                            {
                                connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                                    // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                                    update_property_tree_widget_name( treeWidgetItem, aText );
                                } );
                            }
                            else if ( tRoot == (uint)Module_Type::FEM && tChild == (uint)FEM_Submodule::PHASES && mQTreeWidgetSubChildren[ tRoot ][ tChild ].size() > 0 )
                            {
                                connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                                    // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                                    update_phase_tree_widget_name( treeWidgetItem, aText );
                                } );
                            }
                            else
                            {

                                connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                                    // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                                    update_tree_widget_name( treeWidgetItem, aText );
                                } );
                            }
                            break;
                        }
                    }
                }
            }
            else
            {
                QList< uint > tCurrentIndex = tCurrentWidget->getIndex();
                mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).add_parameter_list( create_parameter_list( (Module_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], 0 ) );
                tCurrentWidget->add_elements( mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] )( mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).size() - 1 ) );
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
                if ( tCurrentWidget->getSubFormType() != -1 )
                {
                    mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->removeSubFormCountProps( tCurrentWidget->getSubFormType() );
                }
                else
                {
                    mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->removeSubFormCountProps();
                }

                // If the first form is selected then the user is redirected to the main form
                if ( tCurrentIndex[ 2 ] == 0 && mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).size() == 1 )
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
                mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).erase( tCurrentIndex[ 2 ] );
            }
            else
            {
                QMessageBox::warning( this, "Warning", "Cannot remove parameters here." );
            }
        }
    }

    void Moris_Gui::add_project( uint aRoot, uint aChild, uint aSubChild )
    {
        /**
         * @brief Function to add a project to the tree widget
         * @param uint aRoot, uint aChild, uint aSubChild, const std::string &aName
         * @return NONE
         * @note This function adds a project to the tree widget
         */
        QString tName = QString::fromStdString( get_inner_sub_parameter_list_name( (Module_Type)aRoot, aChild ) ) + " " + QString::number( aSubChild );

        for ( auto &iFindName : mLibrary.get_parameter_lists()( aRoot )(aChild)( aSubChild ) )
        {
            if ( endswith( iFindName.get_name(), "name" ) )
            {
                // Check if the name is empty or whitespace then break
                if ( iFindName.get_parameter().get_value< std::string >().empty() )
                {
                    break;
                }
                tName = QString::fromStdString( iFindName.get_parameter().get_value< std::string >() );
                break;
            }
        }

        if ( mTreeWidgetItems.size() >= (uint)Module_Type::FEM )
        {
            if ( mTreeWidgetChildren[ (uint)Module_Type::FEM ].size() >= (uint)FEM_Submodule::CONSTITUTIVE_MODELS )
            {
                if ( mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size() > 0 )
                {
                    for ( uint i = 0; i < mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size(); i++ )
                    {
                        if ( !mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ][ i ]->isPropertyListSet() )
                        {
                            mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ][ i ]->setPropertyNameList( mPropertyNameList );
                        }
                        // pass the updated property name list to the Moris_Tree_Widget_Item by reference
                    }
                }
            }
        }

        mTreeWidgetChildren[ aRoot ][ aChild ]->setSubFormCheck( true );
        mQTreeWidgetSubChildren[ aRoot ][ aChild ].append( new QTreeWidgetItem( mQTreeWidgetChildren[ aRoot ][ aChild ] ) );
        mQTreeWidgetChildren[ aRoot ][ aChild ]->addChild( mQTreeWidgetSubChildren[ aRoot ][ aChild ].last() );
        mQTreeWidgetSubChildren[ aRoot ][ aChild ].last()->setText( 0, tName );
        mTreeWidgetSubChildren[ aRoot ][ aChild ].append( new Moris_Tree_Widget_Item() );
        mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->setupScrollArea();
        mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->setIndex( { aRoot, aChild, (uint)mQTreeWidgetSubChildren[ aRoot ][ aChild ].size() - 1 } );
        mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->setSubFormType( -1 );
        mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->add_elements( mLibrary.get_parameter_lists()( aRoot )(aChild)( aSubChild ) );
        mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->set_form_visible( false );
        mTreeWidget->setItemWidget( mQTreeWidgetSubChildren[ aRoot ][ aChild ].last(), 0, mTreeWidgetSubChildren[ aRoot ][ aChild ].last() );
        mLayout->addWidget( mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->getScrollArea() );

        if ( mQTreeWidgetChildren[ aRoot ][ aChild ]->text( 0 ).toStdString() == get_submodule_names( Module_Type::FEM )( 0 ) && mQTreeWidgetSubChildren[ aRoot ][ aChild ].size() > 0 )
        {
            mPropertyNameList.resize( mQTreeWidgetSubChildren[ aRoot ][ aChild ].size() );
            for ( uint i = 0; i < mQTreeWidgetSubChildren[ aRoot ][ aChild ].size(); i++ )
            {
                mPropertyNameList[ i ] = mQTreeWidgetSubChildren[ aRoot ][ aChild ][ i ]->text( 0 );
            }
        }

        // Get mWidget from the Moris_Tree_Widget_Item by reference
        QList< QWidget * > &tWidget = mTreeWidgetSubChildren[ aRoot ][ aChild ].last()->mWidget;

        for ( uint it = 0; it < tWidget.size(); it++ )
        {
            if ( endswith( tWidget[ it ]->objectName().toStdString(), "name" ) )
            {
                auto &tLineEdit      = dynamic_cast< Moris_Line_Edit      &>( *tWidget[ it ] );
                auto &treeWidgetItem = mTreeWidgetSubChildren[ aRoot ][ aChild ].last();

                if ( aRoot == (uint)Module_Type::FEM && aChild == (uint)FEM_Submodule::PROPERTIES && mQTreeWidgetSubChildren[ aRoot ][ aChild ].size() > 0 )
                {
                    connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                        // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                        update_property_tree_widget_name( treeWidgetItem, aText );
                    } );
                }
                else if ( aRoot == (uint)Module_Type::FEM && aChild == (uint)FEM_Submodule::PHASES && mQTreeWidgetSubChildren[ aRoot ][ aChild ].size() > 0 )
                {
                    connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                        // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                        update_phase_tree_widget_name( treeWidgetItem, aText );
                    } );
                }
                else
                {

                    connect( &tLineEdit, &QLineEdit::textChanged, this, [ this, treeWidgetItem ]( const QString &aText ) {
                        // Call the rename_tree_widget_item function, passing both the QTreeWidgetItem and the text
                        update_tree_widget_name( treeWidgetItem, aText );
                    } );
                }
                break;
            }
        }
    }

    void Moris_Gui::write_to_xml()
    {
        /**
         * @brief Function to save the current parameters to an XML file
         * @note This function is called when the m_save_to_xml_Button is clicked. It writes the parameters to an XML file.
         */
        QString tFilePath = get_moris_file_path_for_writing();
        mLibrary.print_parameter_receipt( tFilePath.toStdString() );
        QCoreApplication::quit();
    }

    void Moris_Gui::update_tree_widget_name( Moris_Tree_Widget_Item *aItem, const QString &aText )
    {
        /**
         * @brief Function to update the name of the tree widget item
         * @param QTreeWidgetItem *aItem, const QString &aText
         * @return NONE
         * @note This function updates the name of the tree widget item
         */

        // get the index for the Moris_Tree_Widget_Item
        QList< uint > tIndex = aItem->getIndex();
        // set the text of the QTreeWidgetItem
        mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ tIndex[ 2 ] ]->setText( 0, aText );
    }

    void Moris_Gui::update_property_tree_widget_name( Moris_Tree_Widget_Item *aItem, const QString &aText )
    {
        /**
         * @brief Function to update the name of the tree widget item
         * @param QTreeWidgetItem *aItem, const QString &aText
         * @return NONE
         * @note This function updates the name of the tree widget item
         */

        // get the index for the Moris_Tree_Widget_Item
        QList< uint > tIndex = aItem->getIndex();
        // set the text of the QTreeWidgetItem
        mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ tIndex[ 2 ] ]->setText( 0, aText );

        if ( mQTreeWidgetChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ]->text( 0 ).toStdString() == get_submodule_names( Module_Type::FEM )( 0 ) && mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size() > 0 )
        {
            mPropertyNameList.resize( mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size() );
            for ( uint i = 0; i < mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size(); i++ )
            {
                mPropertyNameList[ i ] = mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ i ]->text( 0 );
            }
        }

        if ( mTreeWidgetItems.size() >= (uint)Module_Type::FEM )
        {
            if ( mTreeWidgetChildren[ (uint)Module_Type::FEM ].size() >= (uint)FEM_Submodule::CONSTITUTIVE_MODELS )
            {
                if ( mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size() > 0 )
                {
                    for ( uint i = 0; i < mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ].size(); i++ )
                    {
                        // pass the updated property name list to the Moris_Tree_Widget_Item by reference
                        mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ (uint)FEM_Submodule::CONSTITUTIVE_MODELS ][ i ]->setPropertyNameList( mPropertyNameList );
                    }
                }
            }
        }
    }

    void Moris_Gui::update_phase_tree_widget_name( Moris_Tree_Widget_Item *aItem, const QString &aText )
    {
        /**
         * @brief Function to update the name of the tree widget item
         * @param QTreeWidgetItem *aItem, const QString &aText
         * @return NONE
         * @note This function updates the name of the tree widget item
         */

        // get the index for the Moris_Tree_Widget_Item
        QList< uint > tIndex = aItem->getIndex();
        // set the text of the QTreeWidgetItem
        mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ tIndex[ 2 ] ]->setText( 0, aText );

        if ( mQTreeWidgetChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ]->text( 0 ).toStdString() == get_submodule_names( Module_Type::FEM )( (uint) FEM_Submodule::PHASES ) && mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size() > 0 )
        {
            mPhaseNameList.resize( mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size() );
            for ( uint i = 0; i < mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ].size(); i++ )
            {
                mPhaseNameList[ i ] = mQTreeWidgetSubChildren[ tIndex[ 0 ] ][ tIndex[ 1 ] ][ i ]->text( 0 );
            }
        }

        // Loop over all mWidgets in FEM and wherever there is a mWidget named "phase_name" set_options_list to mPhaseNameList
        if ( mTreeWidgetItems.size() >= (uint)Module_Type::FEM )
        {
            for ( uint i = 0; i < mTreeWidgetChildren[ (uint)Module_Type::FEM ].size(); i++ )
            {
                if ( i != (uint)FEM_Submodule::PHASES )
                {
                    for ( uint iChild = 0; iChild < mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ i ].size(); iChild++ )
                    {
                        QList< QWidget * > &tWidget = mTreeWidgetSubChildren[ (uint)Module_Type::FEM ][ i ][ iChild ]->mWidget;
                        for ( uint it = 0; it < tWidget.size(); it++ )
                        {
                            if ( tWidget[ it ]->objectName().toStdString() == "phase_name" )
                            {
                                auto &tComboBox = dynamic_cast< Moris_Combo_Box & >( *tWidget[ it ] );
                                tComboBox.set_options_list( mPhaseNameList );
                            }
                        }
                    }
                }
            }
        }
    }


    QList< QStringList > Moris_Gui::convert_parameters_to_QStringList( Parameter_List aList )
    {

        /**
         * @brief Function to convert the parameter list to a QList of QStringList
         * @param Parameter_List aParameterList
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
            tQtList[ 0 ].append( QString::fromStdString( it.get_name() ) );
            tQtList[ 1 ].append( QString::fromStdString( it.get_parameter().get_string() ) );
        }

        return tQtList;
    }


}    // namespace moris
#include "main_gui.moc"
