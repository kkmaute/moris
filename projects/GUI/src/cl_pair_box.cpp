#include "cl_pair_box.hpp"

namespace moris
{
    // Constructor for Moris_Pair_Box.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_param: Reference to a moris::Parameter object to be linked with this widget.
    // - a_options: QStringList containing options for the combo box.
    Moris_Pair_Box::Moris_Pair_Box( QWidget *a_parent, moris::Parameter &a_param, const QStringList &a_options )
            : QWidget( a_parent )
            , moris_pair_combo_box( new QComboBox( this ) )
            , moris_pair_line_edit( new QLineEdit( this ) )
            , parameter( a_param )
    {
        // Set up the combo box with provided options
        moris_pair_combo_box->addItems( a_options );
        connect( moris_pair_combo_box, &QComboBox::currentTextChanged, this, &Moris_Pair_Box::on_combo_box_text_changed );

        // Set up the line edit
        connect( moris_pair_line_edit, &QLineEdit::textChanged, this, &Moris_Pair_Box::on_line_edit_text_changed );

        // Layout the widgets in a horizontal layout
        auto *layout = new QHBoxLayout( this );
        layout->addWidget( moris_pair_combo_box );
        layout->addWidget( moris_pair_line_edit );
        setLayout( layout );
    }

    // Slot to handle changes in the combo box text
    // Inputs:
    // - a_text: New text selected in the combo box.
    void Moris_Pair_Box::on_combo_box_text_changed( const QString &a_text )
    {
        emit combo_box_text_changed( objectName(), a_text );
    }

    // Slot to handle changes in the line edit text
    // Inputs:
    // - a_text: New text entered in the line edit.
    void Moris_Pair_Box::on_line_edit_text_changed( const QString &a_text )
    {
        emit line_edit_text_changed( objectName(), a_text );
    }

    // Set the associated moris::Parameter with the combined value of the combo box and line edit
    void Moris_Pair_Box::set_parameter()
    {
        std::pair< std::string, std::string > pair_value;
        pair_value.first  = moris_pair_combo_box->currentText().toStdString();
        pair_value.second = moris_pair_line_edit->text().toStdString();
        parameter.set_value( objectName().toStdString(), pair_value, false );
    }
}    // namespace moris
