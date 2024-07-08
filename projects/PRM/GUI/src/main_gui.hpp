/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main_gui.hpp
 *
 */

#ifndef MAIN_GUI_H
#define MAIN_GUI_H

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

#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_PRM_MIG_Parameters.hpp"
#include "fn_PRM_STK_Parameters.hpp"
#include "fn_PRM_WRK_Parameters.hpp"


class Moris_Gui : public QWidget
{
    Q_OBJECT

  public:
    // Add argument that reads a parameter list
    Moris_Gui( QWidget *parent = nullptr );

  private slots:

    // Explain the inputs and outputs of the functions

    QList< QStringList > get_parameter_list( moris::Parameter_List_Type, int, int );
    QStringList          get_outer_sub_parameter_list( moris::Parameter_List_Type aModule );

    void set_form_visible( int, int, int, bool );
    void add_elements( int, int, int );

    // CM_Struc_Linear_Isotropic::CM_Struc_Linear_Isotropic();

    void setup_scroll_widget( int, int, int );
    void add_more_props();
    void remove_props();

    QList< int > get_tree_index( QTreeWidgetItem * );

    void parameter_selected( QTreeWidgetItem *, QTreeWidgetItem * );

    QList< QStringList > convert_parameters_to_QStringList( moris::Parameter_List );

  private:
    // Layout objects
    /*
    mLayout (QHBoxLayout) -> This is the main layout in the GUI
    mSidePanel (QVBoxLayout) -> This is the layout on the leftmost side of mLayout to show the TreeWidget for navigation
    mFormLayout -> This is a 2D list that navigates to the corresponding submodules
      For e.g. mFormLayout[root][child][subchild]
                root -> corresponds to the Project
                child -> corresponds to the submodule in a certain peoject (root)
                subchild -> some submodules (child) may have inner sub modules and this corresponds to that
    mAddButton (QPushButton) -> Add button to add more parameters in a certain form
    mRemoveButton (QPushButton) -> Remove button to remove parameters in a certain form
    */
    QHBoxLayout                             *mLayout    = new QHBoxLayout( this );
    QVBoxLayout                             *mSidePanel = new QVBoxLayout();
    QList< QList< QList< QFormLayout * > > > mFormLayout;
    QPushButton                             *mAddButton    = new QPushButton;
    QPushButton                             *mRemoveButton = new QPushButton;

    // Scrolling objects

    /*
    Adding the mFormLayout object to its corresponding mScrollWidget and mScrollArea enables scrolling for the form
    */

    QList< QList< QList< QWidget * > > >     mScrollWidget;
    QList< QList< QList< QScrollArea * > > > mScrollArea;

    QStringList mProjectNames;

    // Tree Widget items for navigation

    /*
    mTreeWidget (QTreeWidget) -> Tree widget in mSidePanel to enable navigation through the GUI
    mTreeWidgetItems -> This list of QTreeWidgetItems (QTreeWidgetItems are items that are added to a QTreeWidget (mTreeWidget),
      this keeps track of the root (Project) added to the treeWidget
    mTreeWidgetChidren -> This 2D list of QTreeWidgetItems that keep track of the children (sub-module) items of a root (Project) added to mTreeWidget
    mTreeWidgetSubChildren -> This is a 3D List of QTreeWidgetItems that keeps track of the SubChildren (inner sub-modules) added to mTreeWidget

    mOldSelection -> QList of 3 ints keeps track of the oldSelection on the chance that the user selects
      a QTreeWidgetItem that does not have a corresponding QFormLayout
    */
    QTreeWidget                                 *mTreeWidget = new QTreeWidget();
    QList< QTreeWidgetItem * >                   mTreeWidgetItems;
    QList< QList< QTreeWidgetItem * > >          mTreeWidgetChildren;
    QList< QList< QList< QTreeWidgetItem * > > > mTreeWidgetSubChildren;
    QList< int >                                 mOldSelection;

    // Element related objects

    /*
    mLineEdit -> This is a 4D QList of QLineEdit objects, the first 3 indices correspond to the root, children, subchildren
      and the last index keeps track of the lineEdit fields inside the following QFormLayout
    mCountProps -> This is a 3D QList of ints that keeps track of the number of parameters in a form
    mOPTAlgorithmComboBox (QComboBox) -> this combo box keeps all the possible algorithms that can be selected.
      OPT (Project) -> Algorithm (child) has a different structure,
      normally the mAddButton, adds parameters (A bunch of rows containing QLineEdit objects) in the same form but in this case,
      an algorithm is selected from the combo box and then in the Algorithms sub-modules, there is a new inner sub-module is added based on the selection
    mOPTAlgorithmComboBox (QComboBox) -> this combo box keeps all the possible algorithms that can be selected.
      GEN (Project) -> Geometries (child) has the same functionality as OPT (Project) -> Algorithm (child)
    */
    QList< QList< QList< QList< QLineEdit * > > > > mLineEdit;
    QList< QList< QList< int > > >                  mCountProps;

    
    moris::Parameter_List mParam = moris::prm::create_gcmma_parameter_list();
    QComboBox                                      *mOPTAlgorithmComboBox  = new QComboBox;
    QComboBox                                      *mGENGeometriesComboBox = new QComboBox;
    QList< QList< QComboBox * > >                   mComboBox;
};
#endif
