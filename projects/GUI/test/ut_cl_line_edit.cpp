// MGUIT431: Moris Gui Test Implement
#include "ut_cl_line_edit.hpp"


#include <QSignalSpy>
#include "../src/cl_line_edit.hpp"

using namespace moris;

void ut_cl_line_edit::initTestCase()
    {
        QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
    }
void ut_cl_line_edit::testTextChanged()
    {
        // Set up a base parameter with no external validators to test on
        Parameter testParameter(std::string("initial"),
                                Entry_Type::FREE, // no need for external parameters
                                "", //not necessary because no external parameters
                                Module_Type::END_ENUM, 
                                0);
        // create a line edit with no parent and our testParameter
        moris::Moris_Line_Edit testLineEdit(nullptr, testParameter);
        testLineEdit.setObjectName("myEdit"); // a default call of a QObject
        
        QSignalSpy textChangedSpy(&testLineEdit, &moris::Moris_Line_Edit::textChanged);
        QVERIFY(textChangedSpy.isValid());

        testLineEdit.clear();
        QTest::keyClicks(&testLineEdit, "textchangeconfirmed");

        QVERIFY(textChangedSpy.count() >= 1); //check that at least one signal has occured
        
        auto textChangeSignal = textChangedSpy.takeLast();
        QCOMPARE(textChangeSignal.at(0).toString(), QString("myEdit"));
        QCOMPARE(textChangeSignal.at(1).toString(), QString("textchangeconfirmed"));
        QCOMPARE(testParameter.get_value<std::string>(), std::string("textchangeconfirmed"));
    }




