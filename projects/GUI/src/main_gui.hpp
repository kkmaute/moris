/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main_gui.hpp
 *
 */

#pragma once

#include <QApplication>
#include <QWidget>
#include <QFormLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QComboBox>
#include <QScrollArea>
#include <QTreeWidget>
#include <QList>

#include <iostream>

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src

#include "cl_Library_IO.hpp"    // MRS/IOS/src
#include "moris_tree_widget_item.hpp"

// Using moris namespace

namespace moris
{
    class Moris_Gui : public QWidget
    {
        Q_OBJECT

        //--------------------------------------------------------------------------------------------------------------

        /** Layout objects
         * mLayout (QHBoxLayout) -> This is the main layout in the GUI
         * mSidePanel (QVBoxLayout) -> This is the layout on the leftmost side of mLayout to show the TreeWidget for navigation
         * mAddButton (QPushButton) -> Add button to add more parameters in a certain form
         * mRemoveButton (QPushButton) -> Remove button to remove parameters in a certain form
         */

        QHBoxLayout *mLayout       = new QHBoxLayout( this );
        QVBoxLayout *mSidePanel    = new QVBoxLayout();
        QPushButton *mAddButton    = new QPushButton;
        QPushButton *mRemoveButton = new QPushButton;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * mProjectNames (QStringList) -> This QStringList keeps the names of the projects
         */
        QStringList mProjectNames;

        //--------------------------------------------------------------------------------------------------------------

        /** Element related objects
         * mTreeWidgetItems (QList< Moris_Tree_Widget_Item * >) -> This list of Moris_Tree_Widget_Items (Moris_Tree_Widget_Item is a class that is used to create a form in the GUI)
         * mTreeWidgetChildren (QList< QList< Moris_Tree_Widget_Item * > >) -> This 2D list of Moris_Tree_Widget_Items that keep track of the children (sub-module) items of a root (Project) added to mTreeWidget
         * mTreeWidgetSubChildren (QList< QList< QList< Moris_Tree_Widget_Item * > > >) -> This is a 3D List of Moris_Tree_Widget_Items that keeps track of the sub-forms (inner sub-modules) added to mTreeWidget
         */
        QList< Moris_Tree_Widget_Item * >                   mTreeWidgetItems;
        QList< QList< Moris_Tree_Widget_Item * > >          mTreeWidgetChildren;
        QList< QList< QList< Moris_Tree_Widget_Item * > > > mTreeWidgetSubChildren;

        //--------------------------------------------------------------------------------------------------------------

        /** Tree Widget items for navigation
         * mTreeWidget (QTreeWidget) -> Tree widget in mSidePanel to enable navigation through the GUI
         * mQTreeWidgetItems -> This list of QTreeWidgetItems (QTreeWidgetItems are items that are added to a QTreeWidget (mTreeWidget),
         *   this keeps track of the root (Project) added to the treeWidget
         * mQTreeWidgetChidren -> This 2D list of QTreeWidgetItems that keep track of the children (sub-module) items of a root (Project) added to mTreeWidget
         * mQTreeWidgetSubChildren -> This is a 3D List of QTreeWidgetItems that keeps track of the sub-forms (inner sub-modules) added to mTreeWidget
         * mOldItem -> This is a pointer to the old item selected in the treeWidget so that if a Moris_Tree_Widget_Item is selected with no form associated, this can keep track of the old item
         */
        QTreeWidget                                 *mTreeWidget = new QTreeWidget();
        QList< QTreeWidgetItem * >                   mQTreeWidgetItems;
        QList< QList< QTreeWidgetItem * > >          mQTreeWidgetChildren;
        QList< QList< QList< QTreeWidgetItem * > > > mQTreeWidgetSubChildren;

        QTreeWidgetItem *mOldItem;

        Vector< Vector< Vector< Parameter_List > > > mParameterLists;

      public:
        // Add argument that reads a parameter list

        /**
         * @brief Constructor for the Moris_Gui class
         * @param QWidget *parent
         * @note This constructor initializes the GUI, sets up the layout, the layout and elements related objects and connects the signals and slots
         */
        Moris_Gui( QWidget *parent = nullptr );

      private slots:

        /**
         * @brief Handles enabling/disabling of Moris_Tree_Widget_Items when a QTreeWidgetItem is selected
         * @note This function is called when a QTreeWidgetItem is selected in the mTreeWidget (It is the slot for itemSelectionChanged() signal of mTreeWidget)
         */
        void parameter_selected();

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief Function to add more parameters when the mAddButton is clicked
         * @note This function is called when the mAddButton is clicked. It adds more parameters in a form or adds a new sub-form if the Moris_Tree_Widget_Item has a sub-form.
         */
        void add_more_props();

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief Function to remove parameters when the mRemoveButton is clicked
         * @note This function is called when the mRemoveButton is clicked. It removes parameters in a form or removes a sub-form if the Moris_Tree_Widget_Item has sub-forms.
         */
        void remove_props();

        //--------------------------------------------------------------------------------------------------------------

      private:
        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief Function to convert the parameter list to a QList of QStringList
         * @param moris::Parameter_List aParameterList
         * @return QList< QStringList >
         * @note This function converts the MORIS parameter list to a QList of QStringList where the 0th index of the list gives the key of the parameter_list and the 1st index gives the value of the parameter_list
         */

        QList< QStringList > convert_parameters_to_QStringList( moris::Parameter_List );
    };


    Vector< Vector< Vector< Parameter_List > > > read();
    /**
     * @brief Function to get the parameter list for a given module, child, and sub-child
     * @param moris::Parameter_List_Type aModule -> root index
     * @param uint aChild -> child index
     * @param uint aSubChild -> sub-child index (if applicable)
     * @return Parameter_List
     * @note This function returns a List of QStringList where the 0th index of the list gives the key of the parameter_list and the 1st index gives the value of the parameter_list
     */
    

}    // namespace moris