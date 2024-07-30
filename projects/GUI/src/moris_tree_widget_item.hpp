#ifndef MORIS_TREE_WIDGET_ITEM_H
#define MORIS_TREE_WIDGET_ITEM_H

#include <QTreeWidgetItem>
#include <QFormLayout>
#include <QScrollArea>
#include <QWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QMessageBox>

#include "moris_line_edit.hpp"
#include "moris_combo_box.hpp"
#include "moris_bool_combo_box.hpp"
#include "moris_double_spin_box.hpp"
#include "moris_int_spin_box.hpp"

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
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"

class Moris_Tree_Widget_Item : public QWidget
{
    Q_OBJECT

  public:
    /**
     * @brief Construct a new Moris_Tree_Widget_Item object that is inherited from QWidget
     * @param QWidget* parent
     * @note trivial constructor
     */
    Moris_Tree_Widget_Item( QWidget *parent = nullptr );
    ~Moris_Tree_Widget_Item();

    /**
     * @brief Function to set up the scroll functionality for any form layout
     * @param NONE
     * @return NONE
     */

    void setupScrollArea();

    /**
     * @brief Function to get the associated scroll area with the Moris_Tree_Widget_Item (to add to the main layout)
     * @param NONE
     * @return QScrollArea*
     */
    QScrollArea *getScrollArea();

    /**
     * @brief Function to set the items passed in to the combo box
     * @param QStringList aItems
     * @return NONE
     */
    void setComboBoxItems( QStringList aItems );

    /**
     * @brief Function to get the associated combo box with the Moris_Tree_Widget_Item
     * @param NONE
     * @return QComboBox*
     */
    QComboBox *getComboBox();

    /**
     * @brief Function to set the sub form check (if the form has a sub form)
     * @param bool aSubFormCheck
     * @return NONE
     */

    void setSubFormCheck( bool aSubFormCheck );
    /**
     * @brief Function to check if the form has a sub form
     * @param NONE
     * @return bool
     * @retval true if this form has a sub form
     */

    bool hasSubForm();

    /**
     * @brief Function to check if the form is a sub form
     * @param NONE
     * @return bool
     * @retval true if the form is a sub form
     */
    bool isSubForm();

    /**
     * @brief Function to check if the Moris_Tree_Widget_Item has a form (for e.g. Projects such as FEM do not have an associated form)
     * @param NONE
     * @return bool
     * @retval true if the Moris_Tree_Widget_Item has a form
     */

    bool hasForm();

    // Keep special forms on hold
    /**
     * @brief Function to set the special form status
     * @param bool aSpecialForm
     * @return NONE
     */
    void setSpecialFormStatus( bool aSpecialForm );

    /**
     * @brief Function to check if the form is one of the special forms
     * @param NONE
     * @return bool
     */
    bool isSpecialForm();

    /**
     * @brief Function to set the count of properties in the form
     * @param int aCountProps
     * @return NONE
     */

    void setCountProps( int aCountProps );

    /**
     * @brief Function to get the count of properties in the form
     * @param NONE
     * @return int
     */
    int getCountProps();

    /**
     * @brief Function to set the index of the form (used to navigate mTreeWidgetItems)
     * @param QList < uint > aIndex
     * @return NONE
     * @note This is a temporary function to set the index of the form in the mTreeWidgetItems list. This may not be needed later.
     */
    void setIndex( QList< uint > aIndex );

    /**
     * @brief Function to get the index of the form
     * @param NONE
     * @return QList < uint >
     * @note This is a temporary function to get the index of the form in the mTreeWidgetItems list. This may not be needed later.
     */
    QList< uint > getIndex();

    /**
     * @brief Function to add the number of sub-forms of a certain type.
     * @param NONE
     * @return NONE
     * @note If a form has associated sub forms, it adds the number of sub-forms of the type selected in the mComboBox to keep track of the number of properties of each type.
     */

    void addSubFormCountProps();

    /**
     * @brief Function to remove the count of a sub-form type from mSubFormCountPropsVec
     * @param int aIndex
     * @return NONE
     * @note If a form has associated sub forms, this function will remove the count of properties in the sub form to keep track of the numnber of properties of each type
     */

    void removeSubFormCountProps( int aIndex );

    /**
     * @brief Function to get the number of sub-forms associated with the current selection of mComboBox.
     * @param NONE
     * @return uint
     */
    uint getSubFormCountProps();

    /**
     * @brief Function to set the type of the sub form (from the index of the mComboBox)
     * @param int aSubFormType
     * @return NONE
     * @note If a form has an associated sub form, this function stores the index of the mComboBox selection to keep track of the sub form type.
     */
    void setSubFormType( int aSubFormType );

    /**
     * @brief Function to get the sub form type
     * @param NONE
     * @return int
     * @note If a form has an associated sub form, this function returns the type of the sub form which is the index in the mComboBox
     */
    int getSubFormType();

    /**
     * @brief Function to set the form visible
     * @param bool aVisible
     * @return NONE
     * @note This function sets the visibility of the form and it's associated elements. This is used to hide/show the form based the selection in mTreeWidget.
     */
    void set_form_visible( bool aVisible );

    /**
     * @brief Function to add elements to the form
     * @param QList< QStringList > aParameters
     * @return NONE
     * @note This function adds the parameters passed to it to the form.
     */
    //   void add_elements(QList< QStringList > aParameters);
    void add_elements( moris::Parameter_List &aParameters );

    /**
     * @brief Function to remove elements from the form
     * @param NONE
     * @return NONE
     * @note This function removes one parameter worth of elements and reorganizes mCountProps.
     */
    void remove_elements();

  private:
    //--------------------------------------------------------------------------------------------------------------

    /** Layout and Scrolling Objects
     * mFormLayout -> Associated form layout for Moris_Tree_Widget_Item
     * mScrollArea -> Associated scroll area for Moris_Tree_Widget_Item. This handles the scrolling functionality for the form layout.
     * mScrollWidget -> Associated scroll widget for Moris_Tree_Widget_Item. This is the widget that is scrolled in the mScrollArea and keeps the form layout
     */
    QWidget     *mScrollWidget;
    QFormLayout *mFormLayout;
    QScrollArea *mScrollArea;

    //--------------------------------------------------------------------------------------------------------------

    /** Checkers and Counters
     * mHasSubFormCheck -> Check if the Moris_Tree_Widget_Item has a sub form
     * mIsSubFormCheck -> Check if the Moris_Tree_Widget_Item is a sub form
     * mCheckForm -> Check if Moris_Tree_Widget_Item has an associated form
     * mIsSpecialForm -> Check if the form is a special form
     * mCountProps -> Count of properties in the form
     * mSubFormType -> Type of sub form (keeps the index of mComboBox)
     * mSubFormCountPropsVec -> List keeps the index of the mComboBox and each index keeps the count of sub forms of that type
     */

    bool          mHasSubFormCheck = false;
    bool          mIsSubFormCheck  = false;
    bool          mCheckForm       = false;
    // Keep on hold
    bool          mIsSpecialForm   = false;
    
    int           mCountProps      = 0;
    int           mSubFormType;
    QList< uint > mSubFormCountPropsVec;
    uint         mSubFormCountProps = 0;


    //--------------------------------------------------------------------------------------------------------------

    /**
     * mWidget -> List of widgets in the form layout (including QLineEdits and QComboBoxes for selection items)
     * mComboBox -> Combo box to select sub forms
     *
     */
    QList< QWidget * > mWidget;
    QComboBox         *mComboBox = new QComboBox();

    QList< uint > mIndex;
};

#endif