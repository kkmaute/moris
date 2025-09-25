#include "ut_cl_double_spin_box.hpp"
#include "moris_typedefs.hpp"

#include <QSignalSpy>
#include "../src/cl_double_spin_box.hpp"

using namespace moris;

void ut_cl_double_spin_box::initTestCase()
{
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
}
void ut_cl_double_spin_box::testInitialSigned()
{
    Parameter testParameter(static_cast<real>(-5.5L),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);

    Moris_Double_Spin_Box testSigned(nullptr, testParameter);
    QCOMPARE(testSigned.value(), static_cast<real>(-5.5L));
    QVERIFY(testSigned.minimum() < static_cast<real>(0));
}
void ut_cl_double_spin_box::testInitialUnsigned()
{
    Parameter testParameter(static_cast<real>(8.8L),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);

    Moris_Double_Spin_Box testUnsigned(nullptr, testParameter);
    QCOMPARE(testUnsigned.value(), static_cast<real>(8.8L));
    QVERIFY(testUnsigned.minimum() < static_cast<real>(0));
}
void ut_cl_double_spin_box::testvaluechanged()
{
    Parameter testParameter(static_cast<real>(0),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);
    Moris_Double_Spin_Box testValueChanged(nullptr, testParameter);
    testValueChanged.setObjectName("testDoubleSpinWorking");

    QSignalSpy setValueSpy(&testValueChanged, &Moris_Double_Spin_Box::value_changed);
    QVERIFY(setValueSpy.isValid());

    testValueChanged.setValue(static_cast<real>(15.5L));

    QCOMPARE(setValueSpy.count(), 1);
    auto spyArgs = setValueSpy.takeFirst();

    QCOMPARE(spyArgs.at(0).toString(), QString("testDoubleSpinWorking"));
    QCOMPARE(spyArgs.at(1).toDouble(), static_cast<double>(15.5L));

    QCOMPARE(testParameter.get_value<real>(), static_cast<real>(15.5L));
}

void ut_cl_double_spin_box::testReadOnlyParameter()
{
    Parameter testParameter(static_cast<real>(3.3),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);
    testParameter.set_value("lockedParam",static_cast<real>(3.3), true);


    Moris_Double_Spin_Box testReadOnly(nullptr, testParameter);
    testReadOnly.setReadOnly(true);
}