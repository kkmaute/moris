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
            , moris_pair_line_edit( new QLineEdit( a_parent ) )
            , moris_pair_line_edit_2( new QLineEdit( a_parent ) )
            , mParameter( a_param )
    {
        // Set up the combo box with provided options
        // moris_pair_combo_box->addItems( a_options );
        // connect( moris_pair_combo_box, &QComboBox::currentTextChanged, this, &Moris_Pair_Box::on_combo_box_text_changed );
        connect( moris_pair_line_edit, &QLineEdit::textChanged, this, &Moris_Pair_Box::on_line_edit_text_changed );
        connect( moris_pair_line_edit_2, &QLineEdit::textChanged, this, &Moris_Pair_Box::on_line_edit_text_changed );


        // Set up the line edit
        // connect( moris_pair_line_edit, &QLineEdit::textChanged, this, &Moris_Pair_Box::on_line_edit_text_changed );
        mPairValue = mParameter.get_value< std::pair< std::string, std::string > >();
        moris_pair_line_edit->setText( QString::fromStdString( mPairValue.first ) );
        moris_pair_line_edit_2->setText( QString::fromStdString( mPairValue.second ) );

        // Layout the widgets in a horizontal layout
        auto *layout = new QHBoxLayout( this );
        // layout->addWidget( moris_pair_combo_box );
        layout->addWidget( moris_pair_line_edit );
        layout->addWidget( moris_pair_line_edit_2 );

        setLayout( layout );
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
        mPairValue.first  = moris_pair_line_edit->text().toStdString();
        mPairValue.second = moris_pair_line_edit_2->text().toStdString();
        mParameter.set_value( objectName().toStdString(), mPairValue, false );
    }

    // Set the associated Parameter with the combined value of the combo box and line edit
    // void Moris_Pair_Box::set_parameter()
    // {
    //     std::pair< std::string, std::string > mPairValue;
    //     mPairValue.first  = moris_pair_combo_box->currentText().toStdString();
    //     mPairValue.second = moris_pair_line_edit->text().toStdString();
    //     mParameter.set_value( objectName().toStdString(), mPairValue, false );
    // }
}    // namespace moris
