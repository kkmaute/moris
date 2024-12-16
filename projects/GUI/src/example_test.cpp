#include <QtTest>
#include "main_gui.hpp"

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test.cpp
 *
 */

#include <catch.hpp>
#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"
#include "main_gui.hpp"

using namespace moris;

//---------------------------------------------------------------

extern "C" void
check_gui_results(const std::string &aFileName, uint aTestCaseIndex)
{
    MORIS_LOG_INFO(" ");
    MORIS_LOG_INFO("Checking GUI Results - Test Case %d.", aTestCaseIndex);
    MORIS_LOG_INFO(" ");

    // Open file and validate existence
    REQUIRE(QFile::exists(QString::fromStdString(aFileName)));

    // Add more specific checks as needed
    MORIS_LOG_INFO("File %s exists.", aFileName.c_str());
}

//---------------------------------------------------------------

TEST_CASE("Moris GUI Tests", "[moris],[gui],[example]")
{
    // Initialize QApplication and GUI
    int argc = 0;
    char* argv[] = {};
    QApplication app(argc, argv);
    
    // Test Case: Load and validate GUI
    SECTION("Initialize and Validate GUI") {
        // Initialize QApplication and Moris_Gui
        // Simulate loading parameters from an existing file
        QString filePath = "test_file.xml";
        QVERIFY(QFile::exists(filePath)); // Ensure the file exists    
        gui = new moris::Moris_Gui( nullptr, "test_file.xml" );
        // Ensure file exists
        check_gui_results(fileName, 0);

        // Additional checks can be performed on GUI elements if required
        MORIS_LOG_INFO("GUI initialized and validated successfully.");
    }
}


    // void testNavigateModifyReset() {
    //     QTreeWidget* treeWidget = gui->findChild<QTreeWidget*>("mTreeWidget");
    //     QVERIFY(treeWidget); // Ensure the tree widget exists

    //     QTreeWidgetItemIterator it(treeWidget);
    //     while (*it) {
    //         QTreeWidgetItem* item = *it;
    //         gui->selectTreeWidgetItem(item);

    //         // Modify and reset parameters
    //         QString originalValue = gui->getCurrentParameterValue();
    //         QString newValue = originalValue + "_modified";
    //         gui->setCurrentParameterValue(newValue);
    //         QCOMPARE(gui->getCurrentParameterValue(), newValue);

    //         gui->setCurrentParameterValue(originalValue);
    //         QCOMPARE(gui->getCurrentParameterValue(), originalValue);

    //         ++it;
    //     }
    // }

    // void testAddAndDeleteParameterList() {
    //     QTreeWidget* treeWidget = gui->findChild<QTreeWidget*>("mTreeWidget");
    //     QVERIFY(treeWidget); // Ensure the tree widget exists

    //     QTreeWidgetItemIterator it(treeWidget);
    //     while (*it) {
    //         QTreeWidgetItem* item = *it;
    //         gui->selectTreeWidgetItem(item);

    //         // Add and delete a new parameter list
    //         gui->addParameterList();
    //         QVERIFY(gui->parameterListAddedSuccessfully());

    //         gui->deleteParameterList();
    //         QVERIFY(gui->parameterListDeletedSuccessfully());

    //         ++it;
    //     }
    // }

    // void testSaveToXML() {
    //     QString filePath = "test_file.xml";
    //     QFile file(filePath);
    //     QVERIFY(file.open(QIODevice::ReadOnly | QIODevice::Text)); // Ensure file can be read
    //     QByteArray originalContent = file.readAll();
    //     file.close();

    //     gui->saveToXML(filePath); // Assumes saveToXML is implemented
    //     QVERIFY(QFile::exists(filePath)); // Ensure file is saved

    //     file.open(QIODevice::ReadOnly | QIODevice::Text);
    //     QByteArray newContent = file.readAll();
    //     file.close();

    //     QVERIFY(originalContent == newContent); // Verify content matches
    // }

    void cleanupTestCase() {
        delete gui;
    }
};

QTEST_MAIN(TestMorisGui)
#include "test_moris_gui.moc"
