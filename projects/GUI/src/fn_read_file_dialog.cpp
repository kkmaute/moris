#include "fn_read_file_dialog.hpp"
#include <QFileDialog>
#include <QStandardPaths>

namespace moris
{
    // Function to open a file dialog and get the selected file path
    QString get_moris_file_path()
    {
        // Open a file dialog to select a file
        QString t_file_path = QFileDialog::getOpenFileName(
                nullptr,                                                             // Parent widget
                "Select File",                                                       // Dialog title
                QStandardPaths::writableLocation( QStandardPaths::HomeLocation ),    // Default directory (home directory)
                "XML Files (*.xml*)"                                                 // Filter for file types (only show XML files)
        );

        // Return the selected file path
        return t_file_path;
    }
}    // namespace moris
