#pragma once

#include "cl_Parameter.hpp"
#include "cl_Variant.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

#include <QWidget>
#include <QComboBox>
#include <QLineEdit>
#include <QStringList>
#include <QFormLayout>

namespace moris
{
    class Parameter;
}

namespace moris
{
    // Moris_Group_Box
    // Custom widget that currently contains a pair of QLineEdit fields.
    // It links to a Parameter object and emits signals when the text or selection changes.
    class Moris_Group_Box : public QWidget
    {
        Q_OBJECT

      public:
        // Constructor
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_param: Reference to a Parameter object to be linked with this widget.
        explicit Moris_Group_Box( QWidget *a_parent, Parameter &a_param, QStringList &a_options );

        // Public member variables
        std::map< std::string, QComboBox * > mWidget;
        QFormLayout                         *mFormLayout = new QFormLayout;

        void set_property_list( QStringList &a_options )
        {
            m_options = a_options;
            for ( auto it = mWidget.begin(); it != mWidget.end(); ++it )
            {
                it->second->clear();
                it->second->blockSignals( true );
                it->second->addItems( m_options );
                it->second->blockSignals( false );
            }
            if ( not m_parameter.get_value< std::string >().empty()) {
                std::stringstream ss( m_parameter.get_value< std::string >() );
                std::string       pair;
                while ( std::getline( ss, pair, ';' ) )
                {
                    std::stringstream pairStream( pair );
                    std::string       tKey, tValue;

                    // Split each pair by comma
                    if ( std::getline( pairStream, tValue, ',' ) && std::getline( pairStream, tKey ) )
                    {
                        mWidget[ tKey ]->blockSignals( true );
                        // Set the current index of the combo box
                        mWidget[ tKey ]->setCurrentIndex( m_options.indexOf( QString::fromStdString( tValue ) ) );
                        mWidget[ tKey ]->blockSignals( false );
                    }
                }
            }
        }

      public slots:

        // Slot to handle text changes in the line edit
        // Inputs:
        // - a_text: New text entered in the line edit.
        // void on_line_edit_text_changed( const QString &a_text );

        void on_combo_box_selection_changed( const int a_index );

        void on_property_selection_changed( const int a_index );


        // Slot to set the associated Parameter with the combined value of the combo box and line edit
        // void set_parameter();

      private:
        Parameter   &m_parameter;    // Reference to the Parameter object
        QStringList &m_options;      // List of options for the combo box
    };

}    // namespace moris
