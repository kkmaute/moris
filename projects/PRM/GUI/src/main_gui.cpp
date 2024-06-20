// #include <QApplication>
// #include <QLineEdit>
// #include <QPushButton>
// #include <QVBoxLayout>
// #include <QWidget>
// #include <iostream>

// #include "cl_Communication_Manager.hpp"    // COM/src
// #include "cl_Logger.hpp"                   // MRS/IOS/src

// #include "fn_PRM_FEM_Parameters.hpp"
// #include "fn_PRM_MSI_Parameters.hpp"
// #include "fn_PRM_SOL_Parameters.hpp"
// #include "fn_PRM_VIS_Parameters.hpp"
// #include "fn_PRM_HMR_Parameters.hpp"
// #include "fn_PRM_GEN_Parameters.hpp"
// #include "fn_PRM_XTK_Parameters.hpp"
// #include "fn_PRM_OPT_Parameters.hpp"

// moris::Comm_Manager gMorisComm;
// moris::Logger       gLogger;

// class Moris_Gui : public QWidget
// {
//     Q_OBJECT

//   public:
//     Moris_Gui( QWidget *parent = nullptr )
//             : QWidget( parent )
//     {
//         // Create the QLineEdit and QPushButton
//         mLineEdit            = new QLineEdit( this );
//         lineEdit2 = new QLineEdit(this);
//         QPushButton *button = new QPushButton( "Store Name", this );

//         // Set up the layout
//         QVBoxLayout *layout = new QVBoxLayout( this );
//         layout->addWidget( mLineEdit );
//         layout->addWidget(lineEdit2);
//         layout->addWidget( button );
//         setLayout( layout );

//         // Connect the button's clicked signal to the appropriate slot
//         connect( button, &QPushButton::clicked, this, &Moris_Gui::storeName );
//     }

//   public slots:
//     void storeName()
//     {
//         // Store the name and print it
//         QString name = mLineEdit->text();
//         std::cout << "Name stored: " << name.toStdString() << std::endl;
//     }

//   private:
//     QLineEdit *mLineEdit;
//     QLineEdit *lineEdit2;
// };

#include "main_gui.h"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

Moris_Gui::Moris_Gui( QWidget *parent )
            : QWidget( parent )
{
  mLayout->addLayout(mSidePanel);
  
  mTreeWidget->setColumnCount(1);

  mProjectNames = {"OPT", "HMR", "STK",
                     "XTK", "GEN", "FEM",
                     "SQL", "MSI", "VIS",
                     "MIG", "WRK", "MORISGENERAL"};

  moris::Parameter_List_Type aModule = moris::Parameter_List_Type::OPT;
  QStringList tStringList;
  QStringList tSubChildrenList;

  for (int iRoot = 0; iRoot < mProjectNames.size(); iRoot++) {

    tStringList = get_outer_sub_parameter_list(aModule);

    mTreeWidgetItems.append(new QTreeWidgetItem);
    mTreeWidgetItems[iRoot]->setText(0,mProjectNames[iRoot]);
    
    QList<QTreeWidgetItem*> tList;
    mTreeWidgetChildren.append(tList);

    QList<QList<QTreeWidgetItem*>> tSubList;
    mTreeWidgetSubChildren.append(tSubList);

    QList<QList<QFormLayout*>> tFormLayout;
    mFormLayout.append(tFormLayout);

    QList<QList<QWidget*>> tScrollWidget;
    mScrollWidget.append(tScrollWidget);

    QList<QList<QScrollArea*>> tScrollArea;
    mScrollArea.append(tScrollArea);

    QList<QList<QList<QLineEdit*>>> tLineEditRoot;
    mLineEdit.append(tLineEditRoot);

    QList<QList<int>> tCountPropsRoot;
    mCountProps.append(tCountPropsRoot);

    for (int iChildren = 0; iChildren < tStringList.size();iChildren++) {

      mTreeWidgetChildren[iRoot].append(new QTreeWidgetItem);
      mTreeWidgetChildren[iRoot][iChildren]->setText(0,tStringList[iChildren]);

      QList<QTreeWidgetItem*> tChildList;
      mTreeWidgetSubChildren[iRoot].append(tChildList);

      QList<QFormLayout*> tFormLayoutChild;
      mFormLayout[iRoot].append(tFormLayoutChild);

      QList<QWidget*> tScrollWidgetChild;
      mScrollWidget[iRoot].append(tScrollWidgetChild);

      QList<QScrollArea*> tScrollAreaChild;
      mScrollArea[iRoot].append(tScrollAreaChild);

      QList<QList<QLineEdit*>> tLineEditChild;
      mLineEdit[iRoot].append(tLineEditChild);

      QList<int> tCountPropsChild;
      mCountProps[iRoot].append(tCountPropsChild);
      
      tSubChildrenList = get_inner_sub_parameter_list(aModule, iChildren);


      if (tSubChildrenList.size() == 0) { 

        mTreeWidgetSubChildren[iRoot][iChildren].append(new QTreeWidgetItem);

        mFormLayout[iRoot][iChildren].append(new QFormLayout);
        mScrollWidget[iRoot][iChildren].append(new QWidget);
        mScrollArea[iRoot][iChildren].append(new QScrollArea);

        QList<QLineEdit*> tLineEditSubChild;
        mLineEdit[iRoot][iChildren].append(tLineEditSubChild);

        int tCountPropsSubChild = 0;
        mCountProps[iRoot][iChildren].append(tCountPropsSubChild);

        mScrollWidget[iRoot][iChildren][0]->setLayout(mFormLayout[iRoot][iChildren][0]);
        mScrollArea[iRoot][iChildren][0]->setHorizontalScrollBarPolicy (Qt::ScrollBarAlwaysOff);
        mScrollArea[iRoot][iChildren][0]->setVerticalScrollBarPolicy (Qt::ScrollBarAsNeeded);
        mScrollArea[iRoot][iChildren][0]->setWidgetResizable (true);
        mScrollArea[iRoot][iChildren][0]->setWidget(mScrollWidget[iRoot][iChildren][0]);
        mLayout->addWidget(mScrollArea[iRoot][iChildren][0]);
        

        if (iRoot == 0 && iChildren == 0) {
          set_form_visible(iRoot, iChildren, 0, true); 
          add_elements(0,0,0); 
          mCountProps[0][0][0] = 1;
        }
        else {
          set_form_visible(iRoot, iChildren, 0, false); 
        }

      }
      else {

        for (int iSubChildren = 0; iSubChildren < tSubChildrenList.size(); iSubChildren++) {
          
          mFormLayout[iRoot][iChildren].append(new QFormLayout);
          mScrollWidget[iRoot][iChildren].append(new QWidget);
          mScrollArea[iRoot][iChildren].append(new QScrollArea);

          QList<QLineEdit*> tLineEditChild;
          mLineEdit[iRoot][iChildren].append(tLineEditChild);

          int tCountPropsChild = 0;
          mCountProps[iRoot][iChildren].append(tCountPropsChild);

          mScrollWidget[iRoot][iChildren][iSubChildren]->setLayout(mFormLayout[iRoot][iChildren][iSubChildren]);
          mScrollArea[iRoot][iChildren][iSubChildren]->setHorizontalScrollBarPolicy (Qt::ScrollBarAlwaysOff);
          mScrollArea[iRoot][iChildren][iSubChildren]->setVerticalScrollBarPolicy (Qt::ScrollBarAsNeeded);
          mScrollArea[iRoot][iChildren][iSubChildren]->setWidgetResizable (true);
          mScrollArea[iRoot][iChildren][iSubChildren]->setWidget(mScrollWidget[iRoot][iChildren][iSubChildren]);
          mLayout->addWidget(mScrollArea[iRoot][iChildren][iSubChildren]);
          

          if (iRoot == 0 && iChildren == 0 && iSubChildren == 0) {
            set_form_visible(iRoot, iChildren, iSubChildren, true); 
            add_elements(0,0,0); 
            mCountProps[0][0][0] = 1; 
          }
          else {
            set_form_visible(iRoot, iChildren, iSubChildren, false); 
          }

          mTreeWidgetSubChildren[iRoot][iChildren].append(new QTreeWidgetItem);
          mTreeWidgetSubChildren[iRoot][iChildren][iSubChildren]->setText(0,tSubChildrenList[iSubChildren]);

        }

        mTreeWidgetChildren[iRoot][iChildren]->addChildren(mTreeWidgetSubChildren[iRoot][iChildren]);
      }


    }

    mTreeWidgetItems[iRoot]->addChildren(mTreeWidgetChildren[iRoot]);

    aModule = static_cast<moris::Parameter_List_Type>(static_cast<int>(aModule) + 1);

  }

  mTreeWidget->addTopLevelItems(mTreeWidgetItems);
  mTreeWidget->setHeaderLabel("Projects");
  mSidePanel->addWidget(mTreeWidget);

  mAddButton->setText("Add");
  mSidePanel->addWidget(mAddButton);

  mRemoveButton->setText("Remove");
  mSidePanel->addWidget(mRemoveButton);

  int tOldSelection = 0;
  OldSelection.append(tOldSelection);
  OldSelection.append(tOldSelection);
  OldSelection.append(tOldSelection);

  connect(mTreeWidget,SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
            this, SLOT(parameter_selected(QTreeWidgetItem*,QTreeWidgetItem*)));
  connect(mAddButton,SIGNAL(clicked()),
            this, SLOT(add_more_props()));
  connect(mRemoveButton,SIGNAL(clicked()),
            this, SLOT(remove_props()));
  
}

QStringList Moris_Gui::get_outer_sub_parameter_list(moris::Parameter_List_Type aModule)
{
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
            "WeakForms",                  // 3
            "QuantitiesOfInterest",       // 4
            "ComputationParameters",      // 5
            "Fields",                     // 6
            //"Materials",                  // 7
            "MaterialModels"              // 8
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
        //tNames = { "Remeshing", "Refinement", "Mapping" };
        break;

    default:
        //MORIS_ERROR( false, "Library_Enums::convert_enum_to_string() - Parameter list type enum unknown." );
        break;
    }

    // check validity of the input
    //uint tNumSubParamLists = tNames.size();

    // retrieve the name for the specific sub-parameter list requested
    return tNames;
}

void Moris_Gui::set_form_visible(int aRoot, int aChildren, int aSubChildren, bool aCheck)
{
    //"aRoot" checks the project index (Root of tree)
    //"aChildren" checks the outer parameter selection (child)
    //"aCheck" checks whether it should turn visible or invisible

    for (int i=0;i<mFormLayout[aRoot][aChildren][aSubChildren]->rowCount();i++){
        mFormLayout[aRoot][aChildren][aSubChildren]->setRowVisible(i,aCheck);

    }
    mScrollWidget[aRoot][aChildren][aSubChildren]->setVisible(aCheck);
    mScrollArea[aRoot][aChildren][aSubChildren]->setVisible(aCheck);
}

QList<QStringList> Moris_Gui::get_parameter_list(moris::Parameter_List_Type aModule, int aChild, int aSubChild)
{

  QList<QStringList> tParameterList;

  switch ( aModule )
  {
  case moris::Parameter_List_Type::OPT:
    switch (aChild)
    {
      case 0:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_opt_problem_parameter_list()); 
        break;
      
      case 1:

        switch (aSubChild)
        {
        case 0:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_opt_interface_parameter_list());

          break;

        case 1:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_opt_interface_manager_parameter_list());

          break;
        
        default:
          break;
        }
      
        break;

      case 2:
        switch (aSubChild)
        {
        case 0:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_gcmma_parameter_list());

          break;

        case 1:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_lbfgs_parameter_list());
          break;

        case 2:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_sqp_parameter_list());

          break;

        case 3:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_sweep_parameter_list());

          break;
        default:
          break;
        }

        break;

      default:
        break;
    }
       //Free
      break;

  case moris::Parameter_List_Type::HMR:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_hmr_parameter_list());

      break;

  case moris::Parameter_List_Type::STK:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_stk_parameter_list());
      break;

  case moris::Parameter_List_Type::XTK:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_xtk_parameter_list());
      break;

  case moris::Parameter_List_Type::GEN:
    switch (aChild)
      {
      case 0:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_gen_parameter_list());
        break;

      case 1:

        switch (aSubChild)
        {
        case 0:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_level_set_geometry_parameter_list(moris::gen::Field_Type::CONSTANT));
          break;
        
        case 1:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_surface_mesh_geometry_parameter_list());
          break;

        case 2:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_voxel_geometry_parameter_list());
          break;
        
        default:
          break;
        }
        break;

      case 2:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_gen_property_parameter_list(moris::gen::Field_Type::CONSTANT));
        break;
      
      default:
        break;
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
    switch (aChild)
      {
      case 0:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_property_parameter_list());
        break;

      case 1:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_constitutive_model_parameter_list());
        break;

      case 2:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_stabilization_parameter_parameter_list());
        break;

      case 3:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_IWG_parameter_list());
        break;

      case 4:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_IQI_parameter_list());
        break;

      case 5:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_computation_parameter_list());
        break;
        
      case 6:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_fem_field_parameter_list());
        break;

      case 7:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_material_model_parameter_list());
        break;
      
      default:
        break;
      }


      break;

  case moris::Parameter_List_Type::SOL:
      tParameterList.resize(8);

      /*
        * For tParameterList[0] need to show different menus
        * based on type of solver
        * Options are: Aztec, Amesos, Belos,Petsc, Eigen_Solver
        * Search for create_linear_algorithm_parameter_list
        * Below is only for Amesos
        */

    switch (aChild)
      {
      case 0:

        switch (aSubChild)
        {
        case 0:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_algorithm_parameter_list_aztec());
          break;
        
        case 1:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_algorithm_parameter_list_amesos());
          break;

        case 2:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_algorithm_parameter_list_belos());
          break;
        
        case 3:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_algorithm_parameter_list_petsc());
          break;
        
        case 4:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_eigen_algorithm_parameter_list());
          break;

        case 5:
        //Need to add ML here
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_algorithm_parameter_list_belos());
          break;

        case 6:
          tParameterList = convert_parameters_to_QStringList(moris::prm::create_slepc_algorithm_parameter_list());
          break;

        default:
          break;
        }
        break;
        
        break;

      case 1:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_linear_solver_parameter_list());
        break;

      case 2:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_nonlinear_algorithm_parameter_list());
        break;

      case 3:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_nonlinear_solver_parameter_list());
        break;

      case 4:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_time_solver_algorithm_parameter_list());
        break;

      case 5:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_time_solver_parameter_list());
        break;
        
      case 6:
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_solver_warehouse_parameterlist());
        break;

      case 7:
        //Need to add Preconditioners
        tParameterList = convert_parameters_to_QStringList(moris::prm::create_material_model_parameter_list());
        break;
      
      default:
        break;
      }

      break;

  case moris::Parameter_List_Type::MSI:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_msi_parameter_list());
      break;

  case moris::Parameter_List_Type::VIS:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_vis_parameter_list()); //

      break; 

  case moris::Parameter_List_Type::MIG:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_mig_parameter_list());

      break;

  case moris::Parameter_List_Type::WRK:
      tParameterList = convert_parameters_to_QStringList(moris::prm::create_wrk_parameter_list());
      break;

  case moris::Parameter_List_Type::MORISGENERAL:
      break;

  default:
      //MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
      break;
  }

  return tParameterList;


}

void Moris_Gui::add_elements(int aRoot, int aChild, int aSubChild)
{
    moris::Parameter_List_Type aModule = static_cast<moris::Parameter_List_Type>(aRoot);
    QList<QStringList> tStringList = get_parameter_list(aModule, aChild, aSubChild);

    int tCounter = mLineEdit[aRoot][aChild][aSubChild].size();

    for (int iElements = 0; iElements < tStringList[0].size(); iElements++) {
        mLineEdit[aRoot][aChild][aSubChild].append(new QLineEdit());
        mLineEdit[aRoot][aChild][aSubChild][iElements]->setPlaceholderText(tStringList[1][iElements]);
        mFormLayout[aRoot][aChild][aSubChild]->
            addRow(tStringList[0][iElements],mLineEdit[aRoot][aChild][aSubChild][iElements+tCounter]);
    }

}

void Moris_Gui::parameter_selected(QTreeWidgetItem * tNewItem, QTreeWidgetItem * tOldItem)
{

    // tNewItem is the new selection by the user
    // tOldItem is the old section by the user

    QList<int> tNewIndex = get_tree_index(tNewItem);
    QList<int> tOldIndex = get_tree_index(tOldItem);

    if (!tNewIndex.isEmpty() && mCountProps[tNewIndex[0]][tNewIndex[1]][tNewIndex[2]] == 0) {
        add_elements(tNewIndex[0],tNewIndex[1],tNewIndex[2]);
        mCountProps[tNewIndex[0]][tNewIndex[1]][tNewIndex[2]]++;
    }

    int tNewRoot, tNewChild, tNewSubChild, tOldRoot, tOldChild, tOldSubChild;

    if (!tOldIndex.isEmpty()) {
        tOldRoot = tOldIndex[0];
        tOldChild = tOldIndex[1];
        tOldSubChild = tOldIndex[2];
    }
    else {
        tOldRoot = 0;
        tOldChild = 0;
        tOldSubChild = 0;

        if (!(tOldIndex.isEmpty() && tNewIndex.isEmpty())) {
            set_form_visible(OldSelection[0],OldSelection[1],OldSelection[2],false);
        }

    }


    if (!tNewIndex.isEmpty()) {
        tNewRoot = tNewIndex[0];
        tNewChild = tNewIndex[1];
        tNewSubChild = tNewIndex[2];
        OldSelection = tNewIndex;
    }
    else {
        return;
    }

    set_form_visible(tOldRoot,tOldChild, tOldSubChild, false);
    set_form_visible(tNewRoot,tNewChild, tNewSubChild, true);

}

QList<int> Moris_Gui::get_tree_index(QTreeWidgetItem *tItem)
{
    QList<int> tIndex;

    for (int iRoot = 0; iRoot < mTreeWidgetChildren.size(); iRoot++) {
        for (int iChildren = 0; iChildren < mTreeWidgetChildren[iRoot].size();iChildren++) {
          if (mTreeWidgetSubChildren[iRoot][iChildren].size() == 1) {
            if (tItem == mTreeWidgetChildren[iRoot][iChildren]) {
                tIndex.append(iRoot);
                tIndex.append(iChildren);
                tIndex.append(0);

                return tIndex;
            }
          }
          else {
            for (int iSubChildren = 0; iSubChildren < mTreeWidgetSubChildren[iRoot][iChildren].size();iSubChildren++) {
              if (tItem == mTreeWidgetSubChildren[iRoot][iChildren][iSubChildren]) {
                tIndex.append(iRoot);
                tIndex.append(iChildren);
                tIndex.append(iSubChildren);

                return tIndex;
            }
          }
          }

            // if (tItem == mTreeWidgetItems[iRoot]) {
            //     tIndex.append(iRoot);
            // }
        }
    }

    return tIndex;
}

QList<QStringList> Moris_Gui::convert_parameters_to_QStringList(moris::Parameter_List tList) {
  QList<QStringList> tQtList;
  QStringList tStringList;
  tQtList.append(tStringList);
  tQtList.append(tStringList);

  for ( auto it = tList.begin(); it != tList.end(); ++it )
  {
      //std::cout << it->first << "\n";
      tQtList[0].append(QString::fromStdString(it->first));
      tQtList[1].append(QString::fromStdString(it->second.get_string()));
      
  }

  return tQtList;
}

QStringList Moris_Gui::get_inner_sub_parameter_list(moris::Parameter_List_Type aModule, int aIndex) {
  QStringList tInner;
  switch ( aModule )
  {
  case moris::Parameter_List_Type::OPT:
      switch (aIndex) {
        case 0:
          tInner = {};
        break;

        case 1:
          tInner = {"Interface", "Interface Manager"};
          break;

        case 2:
          tInner = {"gcmma", "lbfgs", "sqp", "sweep"};
        break;

        default:
        break;
      }

      break;

  case moris::Parameter_List_Type::HMR:
    tInner = {}; 
      break;

  case moris::Parameter_List_Type::STK:
      tInner = {}; 
      break;

  case moris::Parameter_List_Type::XTK:
      tInner = {}; 
      break;

  case moris::Parameter_List_Type::GEN:
      switch (aIndex) {
        case 0:
          tInner = {};
        break;

        case 1:
          tInner = {"Level Set", "Surface Mesh", "Voxel"};
          break;

        case 2:
          tInner = {};
        break;

        default:
        break;
      }
      break;

  case moris::Parameter_List_Type::FEM:
      switch (aIndex) {
        case 0:
          tInner = {};
        break;

        case 1:
          tInner = {};
          break;

        case 2:
          tInner = {};
        break;
        case 3:
          tInner = {};
        break;

        case 4:
          tInner = {};
          break;

        case 5:
          tInner = {};
        break;
        case 6:
          tInner = {};
        break;

        case 7:
          tInner = {};
          break;

        default:
        break;
      }
      break;

  case moris::Parameter_List_Type::SOL:

    switch (aIndex) {
        case 0:
          tInner = {"Aztec", "Amesos", "Belos"
                    ,"PETSC","EigenSolver"
                    , "ML", "Slepc_Solver"};
        break;

        case 1:
          tInner = {};
          break;

        case 2:
          tInner = {};
        break;
        case 3:
          tInner = {};
        break;

        case 4:
          tInner = {};
          break;

        case 5:
          tInner = {};
        break;
        case 6:
          tInner = {};
        break;

        case 7:
          tInner = {};
          break;

        default:
        break;
      }

      break;

  case moris::Parameter_List_Type::MSI:
      tInner = {};
      break;

  case moris::Parameter_List_Type::VIS:
      tInner = {};
      break;

  case moris::Parameter_List_Type::MIG:
      tInner = {};
      break;

  case moris::Parameter_List_Type::WRK:
      tInner = {};
      break;

  case moris::Parameter_List_Type::MORISGENERAL:
      tInner = {};
      break;

  default:
      //MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
      break;
  }

  return tInner;
}

void Moris_Gui::add_more_props() {

  QList<int> tIndex = get_tree_index(mTreeWidget->currentItem());
  add_elements(tIndex[0],tIndex[1],tIndex[2]);
  mCountProps[tIndex[0]][tIndex[1]][tIndex[2]]++;

}

void Moris_Gui::remove_props() {
  QList<int> tIndex = get_tree_index(mTreeWidget->currentItem());

  if (mCountProps[tIndex[0]][tIndex[1]][tIndex[2]] <= 1) {
    QMessageBox::warning(this, "title", "Cannot remove more properties");
  }
  else {

    for (int iRemove = 0; 
      iRemove < mLineEdit[tIndex[0]][tIndex[1]][tIndex[2]].size()/mCountProps[tIndex[0]][tIndex[1]][tIndex[2]];
      iRemove++) {
        mFormLayout[tIndex[0]][tIndex[1]][tIndex[2]]->removeRow
        (mFormLayout[tIndex[0]][tIndex[1]][tIndex[2]]->rowCount()-1);
        mLineEdit[tIndex[0]][tIndex[1]][tIndex[2]].removeLast();
      }

    mCountProps[tIndex[0]][tIndex[1]][tIndex[2]]--;
  }
}


#include "main_gui.moc"
