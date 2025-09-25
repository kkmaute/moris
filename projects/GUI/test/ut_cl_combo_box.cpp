#include "ut_cl_combo_box.hpp"

#include <QSignalSpy>
#include "../src/cl_combo_box.hpp" 

using namespace moris; 

void ut_cl_combo_box::initTestCase()
{
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
}
void ut_cl_combo_box::testInitialEnum()
{
    // Build a parameter the enum validator constructor
    Vector<std::string> options = {"gustave", "lune", "maelle"};
    Parameter testParameter(options);

    Moris_Combo_Box testComboBox(nullptr, testParameter);

    QCOMPARE(testComboBox.count(), int(options.size())); // 3 for the options
    for (int instance = 0; instance < int(options.size()); ++instance)
    {
        QCOMPARE(testComboBox.itemText(instance), QString::fromStdString(options(instance)));
    }
    QCOMPARE(testComboBox.currentIndex(), 0); // default to first item

}
void ut_cl_combo_box::testIndexChanged()
{

    Vector<std::string> options = {"gustave", "lune", "maelle"};
    Parameter testParameter(options);

    Moris_Combo_Box testComboBox(nullptr, testParameter);
    testComboBox.setObjectName("myComboBox");

    QSignalSpy indexChangedSpy(&testComboBox, &Moris_Combo_Box::index_changed);
    QVERIFY(indexChangedSpy.isValid());

    // move to maelle
    testComboBox.setCurrentIndex(2);

    QCOMPARE(indexChangedSpy.count(), 1); // check that exactly one signal has occurred

    auto indexChangeSignal = indexChangedSpy.takeFirst();
    QCOMPARE(indexChangeSignal.at(0).toString(), QString("myComboBox"));
    QCOMPARE(indexChangeSignal.at(1).toInt(), 2);
    // this is important because the parameter is a uint, not a string
    QCOMPARE(testParameter.get_value<uint>(), 2); // should be the index of "maelle"
    
    // check that the current text is "maelle"
    QCOMPARE(testComboBox.currentText(), QString("maelle"));
}
void ut_cl_combo_box::testReadOnlyParameter()
{
    // Set up an enum parameter and then lock it with the set_value true
    Vector<std::string> options = {"gustave", "lune", "maelle"};
    Parameter testParameter(options);

    testParameter.set_value("lock", std::string("gustave"), true);

    QStringList optionsList = {"gustave", "lune", "maelle"};

    Moris_Combo_Box testComboBox(nullptr, testParameter, optionsList);
    QVERIFY(!testComboBox.isEnabled()); // should be disabled since the parameter is locked
}

void ut_cl_combo_box::testOverloadedConstructor()
{
    // Build a parameter the enum validator constructor
    Vector<std::string> options = {"gustave", "lune", "maelle"};
    Parameter testParameter(std::string("lune"), options);

    QStringList optionsList = {"gustave", "lune", "maelle"};

    Moris_Combo_Box testComboBox(nullptr, testParameter, optionsList);

    QCOMPARE(testComboBox.count(), optionsList.size() + 1 ); // 3 + 1 for the options

    for (int instance = 0; instance < optionsList.size(); ++instance)
    {
        QCOMPARE(testComboBox.itemText(instance + 1), optionsList[instance]);
    }
}