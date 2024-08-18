#include "cl_tree_widget_item.hpp"

namespace moris
{
    Moris_Tree_Widget_Item::Moris_Tree_Widget_Item( QWidget *parent )
            : QWidget( parent )

    {
    }

    Moris_Tree_Widget_Item::~Moris_Tree_Widget_Item() = default;

    void Moris_Tree_Widget_Item::setupScrollArea()
    {
        /**
         * @brief Function to set up the scroll functionality for any form layout
         * @param NONE
         * @return NONE
         * @note This function sets up the scroll area and scroll widget for the form layout
         */
        mScrollWidget = new QWidget;
        mScrollArea   = new QScrollArea;
        mFormLayout   = new QFormLayout;
        mScrollWidget->setLayout( mFormLayout );
        mScrollArea->setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
        mScrollArea->setVerticalScrollBarPolicy( Qt::ScrollBarAsNeeded );
        mScrollArea->setWidgetResizable( true );
        mScrollArea->setWidget( mScrollWidget );
        // mLayout->addWidget( mScrollArea );

        mCheckForm = true;
    }


    QScrollArea *Moris_Tree_Widget_Item::getScrollArea()
    {
        return mScrollArea;
    }

    void Moris_Tree_Widget_Item::setComboBoxItems( QStringList aItems )
    {
        /**
         * @brief Function to set the items passed in to the combo box
         * @param QStringList aItems
         * @return NONE
         * @note If the associated form has an associated sub-form then the function sets the argument passed into mComboBox, adds mComboBox to the layout and sets up mSubFormCountPropsVec
         */

        if ( mHasSubFormCheck )
        {
            mComboBox->addItems( aItems );
            mFormLayout->addRow( "Please select: ", mComboBox );
            mSubFormCountPropsVec.resize( aItems.size() );
            mSubFormCountPropsVec.fill( 0 );
        }
        else
        {
            QMessageBox::warning( this, "Warning", "ComboBox is not allowed in this form. If you would like a sub-form, please set setSubFormCheck to true." );
        }
    }

    QComboBox *Moris_Tree_Widget_Item::getComboBox()
    {
        return mComboBox;
    }

    void Moris_Tree_Widget_Item::setSubFormCheck( bool aSubFormCheck )
    {
        mHasSubFormCheck = aSubFormCheck;
    }

    bool Moris_Tree_Widget_Item::hasSubForm()
    {
        return mHasSubFormCheck;
    }

    bool Moris_Tree_Widget_Item::isSubForm()
    {
        return mIsSubFormCheck;
    }

    bool Moris_Tree_Widget_Item::hasForm()
    {
        return mCheckForm;
    }
    void Moris_Tree_Widget_Item::setSpecialFormStatus( bool aIsSpecialForm )
    {
        mIsSpecialForm = aIsSpecialForm;
    }

    bool Moris_Tree_Widget_Item::isSpecialForm()
    {
        return mIsSpecialForm;
    }

    void Moris_Tree_Widget_Item::setCountProps( int aCountProps )
    {
        mCountProps = aCountProps;
    }

    int Moris_Tree_Widget_Item::getCountProps()
    {
        return mCountProps;
    }

    void Moris_Tree_Widget_Item::addSubFormCountProps()
    {
        /**
         * @brief Function to increment the count of properties in the sub-form
         * @param NONE
         * @return NONE
         * @note If the associated form has an associated sub-form then the function increments the count of properties in the sub-form
         */
        if ( mHasSubFormCheck && mComboBox->count() > 0 )
        {
            mSubFormCountPropsVec[ mComboBox->currentIndex() ]++;
        }
        else if ( mComboBox->count() == 0 && mHasSubFormCheck )
        {
            mSubFormCountProps++;
            // QMessageBox::warning( this, "Warning", "Please add items to the ComboBox before incrementing the count." );
        }
        else
        {
            QMessageBox::warning( this, "Warning", "SubFormCountProps is not allowed in this form. If you would like a sub-form, please set mSubFormCheck to true." );
        }
    }

    void Moris_Tree_Widget_Item::removeSubFormCountProps( int aIndex )
    {
        /**
         * @brief Function to decrement the count of properties in the sub-form
         * @param int aIndex
         * @return NONE
         * @note If the associated form has an associated sub-form then the function decrements the count of properties in the sub-form
         */

        if ( mHasSubFormCheck && mComboBox->count() > 0 )
        {
            mSubFormCountPropsVec[ aIndex ]--;
        }
        else if ( mComboBox->count() == 0 && mHasSubFormCheck )
        {
            mSubFormCountProps--;
            // QMessageBox::warning( this, "Warning", "Please add items to the ComboBox before decrementing the count." );
        }
        else
        {
            QMessageBox::warning( this, "Warning", "SubFormCountProps is not allowed in this form. If you would like a sub-form, please set mSubFormCheck to true." );
        }
    }

    uint Moris_Tree_Widget_Item::getSubFormCountProps()
    {
        if ( mHasSubFormCheck && mComboBox->count() > 0 )
        {
            return mSubFormCountPropsVec[ mComboBox->currentIndex() ];
        }
        else
        {
            return mSubFormCountProps;
        }
    }

    void Moris_Tree_Widget_Item::setIndex( QList< uint > aIndex )
    {
        mIndex = aIndex;
    }

    QList< uint > Moris_Tree_Widget_Item::getIndex()
    {
        return mIndex;
    }

    void Moris_Tree_Widget_Item::setSubFormType( int aSubFormType )
    {
        mSubFormType    = aSubFormType;
        mIsSubFormCheck = true;
    }

    int Moris_Tree_Widget_Item::getSubFormType()
    {
        return mSubFormType;
    }

    void Moris_Tree_Widget_Item::set_form_visible( bool aVisible )
    {
        /**
         * @brief Function to set the form visible
         * @param bool aVisible
         * @return NONE
         * @note This function sets the form visible or invisible based on the argument passed
         */

        for ( uint i = 0; i < (uint)mFormLayout->rowCount(); i++ )
        {
            mFormLayout->setRowVisible( i, aVisible );
        }
        mScrollWidget->setVisible( aVisible );
        mScrollArea->setVisible( aVisible );
    }

    void Moris_Tree_Widget_Item::add_elements( Parameter_List &aParameters )
    {
        /**
         * @brief Function to add elements to the form
         * @param QList< QStringList > aParameters
         * @return NONE
         * @note This function adds the parameters passed to it to the form.
         */

        if ( mIsSubFormCheck && mCountProps > 0 )
        {
            QMessageBox::warning( this, "Warning", "Cannot add a form here" );
        }
        else
        {
            uint tCounter = mWidget.size();
            uint tIndex   = 0;

            for ( auto &iElements : aParameters )
            {
                if ( iElements.second.get_entry_type() == Entry_Type::SELECTION )
                {
                    Moris_Combo_Box *tComboBox = new Moris_Combo_Box( mScrollWidget, iElements.second );
                    tComboBox->setObjectName( QString::fromStdString( iElements.first ) );
                    // Moris_Combo_Box *tComboBox = new Moris_Combo_Box();
                    mWidget.append( tComboBox );
                    mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                }
                else
                {
                    if ( iElements.second.index() == variant_index< bool >() )
                    {
                        // Make a new boolean combo box
                        Moris_Bool_Combo_Box *tComboBox = new Moris_Bool_Combo_Box( mScrollWidget, iElements.second );
                        tComboBox->setObjectName( QString::fromStdString( iElements.first ) );
                        mWidget.append( tComboBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                    }
                    else if ( iElements.second.index() == variant_index< uint >() )
                    {
                        Moris_Int_Spin_Box *tSpinBox = new Moris_Int_Spin_Box( mScrollWidget, iElements.second );
                        tSpinBox->setObjectName( QString::fromStdString( iElements.first ) );
                        mWidget.append( tSpinBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                    }
                    else if ( iElements.second.index() == variant_index< sint >() )
                    {
                        Moris_Int_Spin_Box *tDoubleSpinBox = new Moris_Int_Spin_Box( mScrollWidget, iElements.second );
                        tDoubleSpinBox->setObjectName( QString::fromStdString( iElements.first ) );
                        mWidget.append( tDoubleSpinBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                    }
                    else if ( iElements.second.index() == variant_index< real >() )
                    {
                        Moris_Double_Spin_Box *tDoubleSpinBox = new Moris_Double_Spin_Box( mScrollWidget, iElements.second );
                        tDoubleSpinBox->setObjectName( QString::fromStdString( iElements.first ) );

                        mWidget.append( tDoubleSpinBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                    }
                    else if ( iElements.second.index() == variant_index< std::pair < std::string, std::string>> ()) {
                        Moris_Pair_Box *tPairBox = new Moris_Pair_Box( mScrollWidget, iElements.second);
                        tPairBox->setObjectName( QString::fromStdString( iElements.first ) );
                        mWidget.append( tPairBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                    }
                    else
                    {
                        if (iElements.first == "properties") {
                            Moris_Group_Box *tGroupBox = new Moris_Group_Box(mScrollWidget, iElements.second);
                            tGroupBox->setObjectName( QString::fromStdString( iElements.first ) );
                        mWidget.append( tGroupBox );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                        }
                        else {

                        Moris_Line_Edit *tLineEdit = new Moris_Line_Edit( mScrollWidget, iElements.second );
                        tLineEdit->setObjectName( QString::fromStdString( iElements.first ) );
                        // tLineEdit->setParameter(iElements.second);
                        mWidget.append( tLineEdit );
                        mFormLayout->addRow( QString::fromStdString( iElements.first ), mWidget[ tIndex + tCounter ] );
                        }
                    }
                }
                tIndex++;
            }
            mCheckForm = true;
            mCountProps++;
        }
    }

    void Moris_Tree_Widget_Item::remove_elements()
    {
        /**
         * @brief Function to remove elements from the form
         * @param NONE
         * @return NONE
         * @note This function removes one parameter worth of elements and reorganizes mCountProps.
         */

        if ( mCountProps <= 1 )
        {
            QMessageBox::warning( this, "Warning", "Cannot remove more properties." );
        }
        else
        {
            int tIterator = mWidget.size() / mCountProps;

            for ( int iRemove = 0; iRemove < tIterator; iRemove++ )
            {
                mFormLayout->removeRow( mFormLayout->rowCount() - 1 );
                mWidget.removeLast();
            }
            mCountProps--;
        }
    }
}    // namespace moris