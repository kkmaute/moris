#include "ut_cl_bool_combo_box.hpp"

#include <QSignalSpy>
#include "../src/cl_bool_combo_box.hpp"

using namespace moris;

void ut_cl_bool_combo_box::initTestCase()
{
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
}
void ut_cl_bool_combo_box::testInitialTrue()
{
    // Build a parameter with boolean values
    Parameter testParameter(true,
                            Entry_Type::FREE, // no need for external parameters
                            "", // not necessary because no external parameters
                            Module_Type::END_ENUM,
                            0);

    Moris_Bool_Combo_Box testBoolComboBox(nullptr, testParameter);

    QCOMPARE(testBoolComboBox.count(), 2); // true and false options

    QCOMPARE(testBoolComboBox.itemText(0), QString("true"));
    QCOMPARE(testBoolComboBox.itemText(1), QString("false"));

    QCOMPARE(testBoolComboBox.currentIndex(), 0); // default to true
}
void ut_cl_bool_combo_box::testInitialFalse()
{
    // Build a parameter with boolean values
    Parameter testParameter(false,
                            Entry_Type::FREE, // no need for external parameters
                            "", // not necessary because no external parameters
                            Module_Type::END_ENUM,
                            0);

    Moris_Bool_Combo_Box testBoolComboBox(nullptr, testParameter);

    QCOMPARE(testBoolComboBox.count(), 2); // true and false options

    QCOMPARE(testBoolComboBox.itemText(0), QString("true"));
    QCOMPARE(testBoolComboBox.itemText(1), QString("false"));

    QCOMPARE(testBoolComboBox.currentIndex(), 1); // default to false
}

void ut_cl_bool_combo_box::testIndexChanged()
{
    // Set up a base parameter with no external validators to test on
    Parameter testParameter(true,
                            Entry_Type::FREE, // no need for external parameters
                            "", // not necessary because no external parameters
                            Module_Type::END_ENUM,
                            0);

    Moris_Bool_Combo_Box testBoolComboBox(nullptr, testParameter);
    testBoolComboBox.setObjectName("myBoolComboBox"); // a default call of

    QSignalSpy indexChangedSpy(&testBoolComboBox, &Moris_Bool_Combo_Box::index_changed);
    QVERIFY(indexChangedSpy.isValid());
    
    testBoolComboBox.setCurrentIndex(1); // change to false
    
    QCOMPARE(indexChangedSpy.count(), 1); // check that at least one signal has occurred
    
    auto indexChangeSignal = indexChangedSpy.takeLast();
    
    QCOMPARE(indexChangeSignal.at(0).toString(), QString("myBoolComboBox"));
    QCOMPARE(indexChangeSignal.at(1).toInt(), 1); // should be the index of "false"
    
    QCOMPARE(testParameter.get_value<bool>(), false); // should be false now
}

void ut_cl_bool_combo_box::testReadOnlyParameter()
{
    // Set up a base parameter with boolean values
    Parameter testParameter(true,
                            Entry_Type::FREE, // no need for external parameters
                            "", // not necessary because no external parameters
                            Module_Type::END_ENUM,
                            0);

    // Lock the parameter
    testParameter.set_value("lock", false, true);

    Moris_Bool_Combo_Box testBoolComboBox(nullptr, testParameter);
    QVERIFY(!testBoolComboBox.isEnabled()); // should be disabled since the parameter is locked
}
