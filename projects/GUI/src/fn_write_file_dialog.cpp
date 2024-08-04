#include <QFileDialog>
#include <QInputDialog>
#include <QStandardPaths>
#include <QDir>

namespace moris
{
    QString get_moris_file_path_for_writing()
    {
        // Open a file dialog to select or input the name of the XML file
        QString t_file_path = QFileDialog::getSaveFileName(
                nullptr,                                                             // Parent widget
                "Select or Input XML File",                                          // Dialog title
                QStandardPaths::writableLocation( QStandardPaths::HomeLocation ),    // Default directory (home directory)
                "XML Files (*.xml);;All Files (*)"                                   // Filter for file types (only show XML files)
        );

        if ( t_file_path.isEmpty() )
        {
            // User canceled the dialog
            return QString();
        }

        // Append the .xml extension if not already present
        if ( !t_file_path.endsWith( ".xml", Qt::CaseInsensitive ) )
        {
            t_file_path += ".xml";
        }

        // Return the file path
        return t_file_path;
    }

}    // namespace moris
