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
#include <QSignalBlocker>

#include <map>
#include <string>
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
        // Constructor for Moris_Group_Box.
        // Initializes the group box widget, links it to a Parameter object,
        // stores the list of selectable property options, creates the form layout,
        // and loads the initial rows from the linked Parameter.
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_param: Reference to a Parameter object to be linked with this widget.
        // - a_options: QStringList containing options for the combo box.
        // Outputs:
        // - None.
        explicit Moris_Group_Box( QWidget *a_parent, Parameter &a_param, const QStringList &a_options );
        ~Moris_Group_Box() override = default; // Virtual to ensure proper cleanup in derived classes
        // Getter for the associated Parameter object.
        // Returns the Parameter currently linked with this group box.
        // This function checks that mParameter is not null before dereferencing it.
        // Inputs:
        // - None.
        // Outputs:
        // - Reference to the linked Parameter object.
        Parameter &get_parameter();
        // Setter for the associated Parameter object.
        // Reassigns this group box to a new Parameter object and refreshes the displayed rows.
        // This is useful when widgets are dynamically created or rebound to different Parameters.
        // Inputs:
        // - a_parameter: Reference to the Parameter object to link with this widget.
        // Outputs:
        // - None.
        void setParameter( Parameter &a_parameter );
        // Refreshes the group box display from the linked Parameter object.
        // Clears the existing rows, checks for a valid Parameter, then rebuilds
        // the rows either from the serialized Parameter value or from default
        // constitutive model properties if no value is stored.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        void refresh_data_parameter();
        // Setter for the available property option list.
        // Updates the QStringList used to populate each internal combo box,
        // then rebuilds the displayed rows using the current Parameter data.
        // Inputs:
        // - a_options: New list of selectable property names/options.
        // Outputs:
        // - None.
        void set_property_list( const QStringList &a_options );

        

      public slots:

        // Slot to handle changes in the constitutive model type combo box.
        // When the constitutive type changes, this function rebuilds the property
        // rows for the selected constitutive model and writes the new serialized
        // group box state back to the linked Parameter.
        // Inputs:
        // - a_index: Index of the selected constitutive model type.
        // Outputs:
        // - None.
        
        void on_combo_box_selection_changed( int a_index );
        // Slot to handle changes in one of the property selection combo boxes.
        // When a row's combo box selection changes, this function serializes
        // the current group box state and writes it back to the linked Parameter.
        // Inputs:
        // - a_index: Newly selected combo box index.
        // Outputs:
        // - None.
        void on_property_selection_changed( int a_index );

    private:
        // Clears all rows from the group box form layout.
        // Removes every row currently displayed in the form layout and clears
        // the internal widget map so the group box can be rebuilt cleanly.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        void clear_rows();
        // Adds one property-selection row to the group box.
        // Creates a combo box, fills it with the available property options,
        // sets its selected value if one is provided, connects it to the update slot,
        // and adds it to the form layout under the given key.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        void add_row( const std::string &a_key, const QString &a_selectedText = QString() );
        // Builds the default property rows for the group box.
        // Creates a default constitutive model, retrieves its property map,
        // and adds one row for each property required by that model.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        void build_default_rows();
        // Builds property rows from a serialized Parameter string.
        // Parses the stored string representation of the group box state and
        // recreates each row with its saved key and selected value.
        // Expected serialized format:
        // - "value1,key1;value2,key2;..."
        // Inputs:
        // - a_serialized: Serialized string containing saved row selections.
        // Outputs:
        // - None.
        void build_rows_from_serialized( const std::string &a_serialized );
        // Serializes the current group box state into a string.
        // Reads each row's selected combo box value and its associated key,
        // then combines all rows into a single string that can be stored
        // back into the linked Parameter.
        // Output format:
        // - "value1,key1;value2,key2;..."
        // Inputs:
        // - None.
        // Outputs:
        // - String representation of the current group box selections.
        std::string serialize_current_state() const;

    private:
        Parameter   *mParameter = nullptr;    // Reference to the Parameter object
        QStringList mOptions;      // List of options for the combo box

    public:
        std::map<std::string, QComboBox*> mWidget; // Map of property keys to their corresponding combo box widgets
        QFormLayout *mFormLayout = nullptr; // Layout to organize the widgets in the group box
    };

}    // namespace moris
