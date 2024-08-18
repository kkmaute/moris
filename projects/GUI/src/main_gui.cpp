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
    Moris_Gui::Moris_Gui( QWidget *parent )
            : QWidget( parent )
    {
        QString tFilePath = get_moris_file_path();

        // load the parameter list from the xml file
        mLibrary.load_parameter_list( tFilePath.toStdString(), File_Type::XML_FILE );

        // Update mParameterLists (sitting in cl_Library_IO) with the parameter lists from the xml file
        mLibrary.load_parameters_from_xml();

        // Get mParameterLists from cl_Library_IO
        // mLibrary.get_parameter_lists() = mLibrary.get_parameter_lists();

        //  mLayout is the main layout of the GUI
        //  mSidePanel is the side panel layout where the tree widget and buttons are placed

        mLayout->addLayout( mSidePanel );
        mTreeWidget->setColumnCount( 1 );

        // Read from mParamterLists instead of moris functions, put these functions in the read() function

        for ( uint iProject = 0; iProject < (uint)( Parameter_List_Type::END_ENUM ); iProject++ )
        {
            mProjectNames.append( QString::fromStdString( convert_parameter_list_enum_to_string( (Parameter_List_Type)( iProject ) ) ) );
        }

        Parameter_List_Type tModule = Parameter_List_Type::OPT;
        QStringList         tStringList;
        QStringList         tSubChildrenList;

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
            for ( uint iList = 0; iList < get_number_of_sub_parameter_lists_in_module( tModule ); iList++ )
            {
                tStringList.append( QString::fromStdString( get_outer_sub_parameter_list_name( tModule, iList ) ) );
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

                if ( iRoot == (uint)( Parameter_List_Type::OPT ) && iChildren == (uint)( OPT_SubModule::OPTIMIZATION_PROBLEMS ) )
                {
                    // Setting the OPT/Optimization Problems form visible upon construction of the GUI
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mLibrary.get_parameter_lists()( iRoot )(iChildren)( 0 ) );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->set_form_visible( true );
                    mOldItem = mQTreeWidgetChildren[ iRoot ][ iChildren ];
                }
                else if ( iRoot == (uint)( Parameter_List_Type::OPT ) && iChildren == (uint)( OPT_SubModule::ALGORITHMS ) )
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
                else if ( iRoot == (uint)( Parameter_List_Type::GEN ) && iChildren == (uint)( GEN_SubModule::GEOMETRIES ) )
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
                else if ( iRoot == (uint)( Parameter_List_Type::SOL ) && iChildren == (uint)( SOL_SubModule::LINEAR_ALGORITHMS ) )
                {
                    // Adding the SOL/Linear Algorithms to its combo box
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                    // Hold
                    mTreeWidgetChildren[ iRoot ][ iChildren ]->setSpecialFormStatus( true );
                    QStringList tSolverList;
                    for ( uint iSolver = 0; iSolver < (uint)sol::SolverType_String::values.size() - 1; iSolver++ )
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
                else
                {
                    if ( mLibrary.get_parameter_lists()( iRoot )( iChildren ).size() > 1 )
                    {
                        mTreeWidgetChildren[ iRoot ][ iChildren ]->setSubFormCheck( true );
                        for ( uint i = 0; i < mLibrary.get_parameter_lists()( iRoot )( iChildren ).size(); i++ )
                        {
                            add_project( iRoot, iChildren, i );
                        }
                    }
                    else
                    {
                        // Adding the parameter_list elements to all the forms and setting the visibility to false
                        mTreeWidgetChildren[ iRoot ][ iChildren ]->add_elements( mLibrary.get_parameter_lists()( iRoot )(iChildren)( 0 ) );
                        mTreeWidgetChildren[ iRoot ][ iChildren ]->setCountProps( 1 );
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

            tModule = (Parameter_List_Type)( (uint)( tModule ) + 1 );

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

                mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].append( new QTreeWidgetItem( mQTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ] ) );
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].append( new Moris_Tree_Widget_Item() );

                mQTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->addChild( mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last() );

                // Reorganizing the count of the sub-forms from the item selection in the combo box
                // All these depend on the item selection in the combo box
                mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->addSubFormCountProps();
                if ( mTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->getComboBox()->count() == 0 )
                {
                    mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setText( 0, mQTreeWidgetChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ]->text( 0 ) + QString::number( tCurrentWidget->getSubFormCountProps() ) );
                    mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setSubFormType( -1 );
                    mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).push_back( create_parameter_list( (Parameter_List_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], 0 ) );
                }
                else
                {
                    mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setText( 0, tCurrentWidget->getComboBox()->currentText() + QString::number( tCurrentWidget->getSubFormCountProps() ) );
                    mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setSubFormType( tCurrentWidget->getComboBox()->currentIndex() );
                    mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).push_back( create_parameter_list( (Parameter_List_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], tCurrentWidget->getComboBox()->currentIndex() ) );
                }

                // Setting up the scroll area, index, saving type for the sub-forms
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setupScrollArea();
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->setIndex( { tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], (uint)mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].size() - 1 } );

                // Append to mLibrary.get_parameter_lists() and pass it into add_elements
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->add_elements( mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).back() );
                mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->set_form_visible( false );

                // Connecting the QTreeWidgetItem to the Moris_Tree_Widget_Item
                mTreeWidget->setItemWidget( mQTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last(), 0, mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last() );

                // Adding the scroll area to the main layout
                mLayout->addWidget( mTreeWidgetSubChildren[ tCurrentIndex[ 0 ] ][ tCurrentIndex[ 1 ] ].last()->getScrollArea() );
            }
            else
            {
                QList< uint > tCurrentIndex = tCurrentWidget->getIndex();
                mLibrary.get_parameter_lists()( tCurrentIndex[ 0 ] )( tCurrentIndex[ 1 ] ).push_back( create_parameter_list( (Parameter_List_Type)tCurrentIndex[ 0 ], tCurrentIndex[ 1 ], 0 ) );
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
                // If the current widget is not a sub-form, then remove the elements from the form in the Moris_Tree_Widget_Item
                mLibrary.get_parameter_lists()( tCurrentWidget->getIndex()[ 0 ] )( tCurrentWidget->getIndex()[ 1 ] ).pop_back();
                tCurrentWidget->remove_elements();
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
        mLibrary.finalize( tFilePath.toStdString() );
        QCoreApplication::quit();
    }

    bool Moris_Gui::endswith( const std::string &aString, const std::string &aEnding )
    {
        /**
         * @brief Function to check if a string ends with a particular substring
         * @param const std::string &aString, const std::string &aEnding
         * @return bool
         * @note This function checks if a string ends with a particular substring
         */
        if ( aString.size() >= aEnding.size() )
        {
            // Compare the end of the string with the suffix
            return aString.compare( aString.size() - aEnding.size(), aEnding.size(), aEnding ) == 0;
        }
        else
        {
            return false;
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
        QString tName = QString::fromStdString( get_inner_sub_parameter_list_name( (Parameter_List_Type)aRoot, aChild ) ) + " " + QString::number( aSubChild );

        for ( auto &iFindName : mLibrary.get_parameter_lists()( aRoot )(aChild)( aSubChild ) )
        {
            if ( endswith( iFindName.first, "_name" ) )
            {
                tName = QString::fromStdString( iFindName.second.get_value< std::string >() );
                break;
            }
        }
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
    }


    Vector< Vector< Vector< Parameter_List > > > read()
    {

        // Create the 3d vector
        Vector< Vector< Vector< Parameter_List > > > tParameterList;
        tParameterList.resize( (uint)( Parameter_List_Type::END_ENUM ) );
        for ( uint iRoot = 0; iRoot < (uint)( Parameter_List_Type::END_ENUM ); iRoot++ )
        {
            tParameterList( iRoot ).resize( get_number_of_sub_parameter_lists_in_module( (Parameter_List_Type)iRoot ) );
            for ( uint iChild = 0; iChild < get_number_of_sub_parameter_lists_in_module( (Parameter_List_Type)iRoot ); iChild++ )
            {
                if ( ( iRoot == (uint)( Parameter_List_Type::OPT ) && iChild == (uint)( OPT_SubModule::ALGORITHMS ) )
                        || ( iRoot == (uint)( Parameter_List_Type::GEN ) && iChild == (uint)( GEN_SubModule::GEOMETRIES ) )
                        || ( iRoot == (uint)( Parameter_List_Type::SOL ) && iChild == (uint)( SOL_SubModule::LINEAR_ALGORITHMS ) ) )
                {
                }
                else
                {
                    tParameterList( iRoot )( iChild ).push_back( create_parameter_list( (Parameter_List_Type)iRoot, iChild, 0 ) );
                }
            }
        }

        // Resize based on the projects/sub-projects
        // Populate based on the parameter_list
        return tParameterList;
    }


}    // namespace moris
#include "main_gui.moc"
