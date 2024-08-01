#include <QFileDialog>
#include <QInputDialog>
#include <QStandardPaths>
#include <QDir>

QString get_moris_file_path_for_writing()
{
    // Open a file dialog to select or input the name of the XML file
    QString filePath = QFileDialog::getSaveFileName(
            nullptr,                                                             // Parent widget
            "Select or Input XML File",                                          // Dialog title
            QStandardPaths::writableLocation( QStandardPaths::HomeLocation ),    // Default directory (home directory)
            "XML Files (*.xml);;All Files (*)"                                   // Filter for file types (only show XML files)
    );

    if ( filePath.isEmpty() )
    {
        // User canceled the dialog
        return QString();
    }

    // Append the .xml extension if not already present
    if ( !filePath.endsWith( ".xml", Qt::CaseInsensitive ) )
    {
        filePath += ".xml";
    }

    // Return the file path
    return filePath;
}
