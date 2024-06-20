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
    Moris_Gui( QWidget *parent = nullptr );

  private slots:
    QList< QStringList > get_parameter_list( moris::Parameter_List_Type, int, int );
    QStringList          get_outer_sub_parameter_list( moris::Parameter_List_Type aModule );
    QStringList          get_inner_sub_parameter_list( moris::Parameter_List_Type aModule, int aIndex );

    void set_form_visible( int, int, int, bool );

    void add_elements( int, int, int );

    void add_more_props();
    void remove_props();

    QList< int > get_tree_index( QTreeWidgetItem * );

    void parameter_selected( QTreeWidgetItem *, QTreeWidgetItem * );

    QList< QStringList > convert_parameters_to_QStringList( moris::Parameter_List );

  private:
    // Layout objects
    QHBoxLayout                             *mLayout    = new QHBoxLayout( this );
    QVBoxLayout                             *mSidePanel = new QVBoxLayout();
    QList< QList< QList< QFormLayout * > > > mFormLayout;
    QPushButton                             *mAddButton    = new QPushButton;
    QPushButton                             *mRemoveButton = new QPushButton;

    // Scrolling objects
    QList< QList< QList< QWidget * > > >     mScrollWidget;
    QList< QList< QList< QScrollArea * > > > mScrollArea;

    QStringList mProjectNames;

    // Tree Widget items for navigation
    QTreeWidget                                 *mTreeWidget = new QTreeWidget();
    QList< QTreeWidgetItem * >                   mTreeWidgetItems;
    QList< QList< QTreeWidgetItem * > >          mTreeWidgetChildren;
    QList< QList< QList< QTreeWidgetItem * > > > mTreeWidgetSubChildren;
    QList< int >                                 mOldSelection;

    // Element related objects
    QList< QList< QList< QList< QLineEdit * > > > > mLineEdit;
    QList< QList< QList< int > > >                  mCountProps;
};
#endif
