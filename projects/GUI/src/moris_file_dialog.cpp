#include "moris_file_dialog.hpp"
#include <QFileDialog>

// Function to open a file dialog and get the selected file path
QString getMorisFilePath() {
    // Open a file dialog to select a file
    QString filePath = QFileDialog::getOpenFileName(
        nullptr,                    // Parent widget
        "Select File",              // Dialog title
        "",                         // Default directory (empty string means the current directory)
        "XML Files (*.xml*)"        // Filter for file types (only show XML files)
    );
    
    // Return the selected file path
    return filePath;
}
