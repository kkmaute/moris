#include "cl_pair_box.hpp"

// Some code is left commented to implement a pair box to pair a combo box and line edit
namespace moris
{
    // Constructor for Moris_Pair_Box.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_param: Reference to a Parameter object to be linked with this widget.
    // - a_options: QStringList containing options for the combo box.
    Moris_Pair_Box::Moris_Pair_Box( QWidget *a_parent, Parameter &a_param )
            : QWidget( a_parent )
            //, moris_pair_combo_box( new QComboBox( this ) )
            , moris_pair_line_edit( new QLineEdit( this ) )
            , moris_pair_line_edit_2( new QLineEdit( this ) )
            , mParameter( &a_param )
    {
        // setup default layout based on HBox,
        auto *tLayout = new QHBoxLayout( this );

        tLayout->addWidget( moris_pair_line_edit );
        tLayout->addWidget( moris_pair_line_edit_2 );

        tLayout->setContentsMargins( 0, 0, 0, 0 );

        setLayout( tLayout );

        //
        connect( moris_pair_line_edit,
                 &QLineEdit::textChanged, 
                 this, 
                 &Moris_Pair_Box::on_line_edit_text_changed );
        connect( moris_pair_line_edit_2, 
                 &QLineEdit::textChanged, 
                 this, 
                 &Moris_Pair_Box::on_line_edit_text_changed );

        refreshDataParameter();

        if ( mParameter && mParameter->is_locked() )
        {
            moris_pair_line_edit->setReadOnly( true );
            moris_pair_line_edit_2->setReadOnly( true );
        }
        else
        {
            moris_pair_line_edit->setReadOnly( false );
            moris_pair_line_edit_2->setReadOnly( false );
        }
    }
    // Destructor, qt handles cleanup
    Moris_Pair_Box::~Moris_Pair_Box() = default;
    // Getter for the associated Parameter object.
    // Returns the Parameter currently linked with this widget.
    // Inputs:
    // - None.
    // Outputs:
    // - Reference to the linked Parameter object.
    Parameter &Moris_Pair_Box::get_parameter()
    {
        MORIS_ERROR(mParameter, "Moris_Pair_Box::get_parameter() called with null mParameter.");
        return *mParameter;
    }

    // Setter for the associated Parameter object.
    // Reassigns this widget to a new Parameter object and refreshes the displayed data.
    // Inputs:
    // - a_parameter: Reference to the Parameter object to link with this widget.
    // Outputs:
    // - None.
    void Moris_Pair_Box::setParameter( Parameter &a_parameter )
    {
        // qDebug() << "[PairBox::setParameter]"
        //      << "widget =" << this
        //      << "name =" << objectName()
        //      << "old mParameter =" << mParameter
        //      << "new mParameter =" << &a_parameter;
        mParameter = &a_parameter;

        refreshDataParameter();

        if(mParameter && mParameter->is_locked() )
        {
            moris_pair_line_edit->setReadOnly( true );
            moris_pair_line_edit_2->setReadOnly( true );
        }
        else
        {
            moris_pair_line_edit->setReadOnly( false );
            moris_pair_line_edit_2->setReadOnly( false );
        }
    }
    // Refreshes the data displayed in the pair box from the linked Parameter object.
    // Function reads the current std::pair<std::string, std::string> value from
    // the Parameter and places the first and second values into the two line edits.
    // Signals are blocked while setting text to prevent unwanted calls to
    // on_line_edit_text_changed during initialization or rebinding.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    void Moris_Pair_Box::refreshDataParameter()
    {
        QSignalBlocker tBlocker1( moris_pair_line_edit );
        QSignalBlocker tBlocker2( moris_pair_line_edit_2 );

        if( !mParameter)    
        {
            moris_pair_line_edit->setText( QString() );
            moris_pair_line_edit_2->setText( QString() );
            return;
        }

        if( mParameter->index() != variant_index< std::pair< std::string, std::string > >() )
        {
            qWarning() << "Moris_Pair_Box::refreshDataParameter() called with Parameter of wrong type. Expected std::pair<std::string, std::string>.";
            // fallback: show the parameter string representation in the first box and nothing in the second box because we dont know how to safely parse the string
            moris_pair_line_edit->setText( QString::fromStdString( mParameter->get_string() ) );
            
            moris_pair_line_edit_2->clear();

            return;
        }

        mPairValue = mParameter->get_value< std::pair< std::string, std::string > >();
        moris_pair_line_edit->setText( QString::fromStdString( mPairValue.first ) );
        moris_pair_line_edit_2->setText( QString::fromStdString( mPairValue.second ) );
    }
// Slot to handle changes in the combo box text
// Inputs:
// - a_text: New text selected in the combo box.
// void Moris_Pair_Box::on_combo_box_text_changed( const QString &a_text )
// {
//     emit combo_box_text_changed( objectName(), a_text );
// }

// Slot to handle changes in the first line edit text of the pair
// Inputs:
// - a_text: New text entered in the line edit.
void Moris_Pair_Box::on_line_edit_text_changed( const QString &a_text )
{
    Q_UNUSED( a_text );
    if ( !mParameter || mParameter->is_locked() )
    {
        return;
    }

    if ( mParameter->index() != variant_index< std::pair< std::string, std::string > >() )
    {
        qWarning() << "Moris_Pair_Box::on_line_edit_text_changed() called with Parameter of wrong type. Expected std::pair<std::string, std::string>.";
        return;
    }


    mPairValue.first  = moris_pair_line_edit->text().toStdString();
    mPairValue.second = moris_pair_line_edit_2->text().toStdString();
    mParameter->set_value( objectName().toStdString(), mPairValue, false );
    emit pair_changed(objectName(),
                      moris_pair_line_edit->text(),
                      moris_pair_line_edit_2->text());

}

}    // namespace moris
