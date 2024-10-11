#include "cl_group_box.hpp"

namespace moris
{
    // Constructor for Moris_Group_Box.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_param: Reference to a Parameter object to be linked with this widget.
    // - a_options: QStringList containing options for the combo box.
    Moris_Group_Box::Moris_Group_Box( QWidget *a_parent, Parameter &a_param )
            : QWidget( a_parent )
            , m_parameter( a_param )
    {
        // Set up the combo box with provided options
        // moris_pair_combo_box->addItems( a_options );
        // connect( moris_pair_combo_box, &QComboBox::currentTextChanged, this, &Moris_Pair_Box::on_combo_box_text_changed );
        setLayout( mFormLayout );
        if ( m_parameter.get_value< std::string >().empty() )
        {
            fem::CM_Factory tCMFactory;

            // Error in cl_MTK_Mesh_Manager when try to create_CM;
            std::shared_ptr< fem::Constitutive_Model > tCM  = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            std::map< std::string, uint >             &tMap = tCM->get_property_map();
            //std::map< std::string, uint > tMap;
            for ( const auto &iMap : tMap )
            {
                const std::string &tKey = iMap.first;

                // Create a new QLineEdit
                QLineEdit *tLineEdit = new QLineEdit();

                // Add the key and the QLineEdit pointer to tRows
                mWidget[ tKey ] = tLineEdit;

                connect( mWidget[ tKey ], &QLineEdit::textChanged, this, &Moris_Group_Box::on_line_edit_text_changed );

                // Add the key and the QLineEdit to the QFormLayout
                mFormLayout->addRow( QString::fromStdString( tKey ), tLineEdit );
            }
        }
        else
        {
            std::stringstream ss( m_parameter.get_value< std::string >() );
            std::string       pair;
            while ( std::getline( ss, pair, ';' ) )
            {
                std::stringstream pairStream( pair );
                std::string       tKey, tValue;

                // Split each pair by comma
                if ( std::getline( pairStream, tValue, ',' ) && std::getline( pairStream, tKey ) )
                {
                    // Create a new QLineEdit
                    QLineEdit *tLineEdit = new QLineEdit();
                    tLineEdit->setText( QString::fromStdString( tValue ) );

                    // Add the key and the QLineEdit pointer to tRows
                    mWidget[ tKey ] = tLineEdit;

                    connect( mWidget[ tKey ], &QLineEdit::textChanged, this, &Moris_Group_Box::on_line_edit_text_changed );

                    // Add the key and the QLineEdit to the QFormLayout
                    mFormLayout->addRow( QString::fromStdString( tKey ), tLineEdit );
                }
            }
        }
    }

    void Moris_Group_Box::on_line_edit_text_changed( const QString &a_text )
    {
        std::string tResult;
        for ( auto it = mWidget.begin(); it != mWidget.end(); ++it )
        {
            // Append the key
            tResult += it->second->text().toStdString() + ",";

            // Append the text from the QLineEdit
            tResult += it->first;

            // If it's not the last element, add a semicolon
            if ( std::next( it ) != mWidget.end() )
            {
                tResult += ";";
            }
        }
        m_parameter.set_value( objectName().toStdString(), tResult, false );
    }

    void Moris_Group_Box::on_combo_box_selection_changed( const int a_index )
    {
        // Clear the form layout
        for ( auto it = mWidget.begin(); it != mWidget.end(); ++it ) {
                mFormLayout->removeRow( mFormLayout->rowCount() - 1 );
        }
        mWidget.clear();
        
        fem::CM_Factory tCMFactory;

            // Error in cl_MTK_Mesh_Manager when try to create_CM;
            // Some constitutive_types are not available in the create_CM function. Throws runtime error.
            std::shared_ptr< fem::Constitutive_Model > tCM  = tCMFactory.create_CM( (fem::Constitutive_Type) a_index );
            std::map< std::string, uint >             &tMap = tCM->get_property_map();
            //std::map< std::string, uint > tMap;
            for ( const auto &iMap : tMap )
            {
                const std::string &tKey = iMap.first;

                // Create a new QLineEdit
                QLineEdit *tLineEdit = new QLineEdit();

                // Add the key and the QLineEdit pointer to tRows
                mWidget[ tKey ] = tLineEdit;

                connect( mWidget[ tKey ], &QLineEdit::textChanged, this, &Moris_Group_Box::on_line_edit_text_changed );

                // Add the key and the QLineEdit to the QFormLayout
                mFormLayout->addRow( QString::fromStdString( tKey ), tLineEdit );
            }
    }

}    // namespace moris
