#include "ut_cl_int_spin_box.hpp"

#include <QSignalSpy>
#include "../src/cl_int_spin_box.hpp"


using namespace moris;

void ut_cl_int_spin_box::initTestCase()
{
    QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
}

void ut_cl_int_spin_box::testInitialSigned()
{
    Parameter testParameter(static_cast<sint>(-5),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);

    Moris_Int_Spin_Box testSigned(nullptr, testParameter);
    QCOMPARE(testSigned.value(), static_cast<int>(-5));
    QVERIFY(testSigned.minimum() < 0);
}

void ut_cl_int_spin_box::testInitialUnsigned()
{
    Parameter testParameter(static_cast<uint>(8),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);

    Moris_Int_Spin_Box testSigned(nullptr, testParameter);
    QCOMPARE(testSigned.value(),8);
    QCOMPARE(testSigned.minimum(), 0);
}

void ut_cl_int_spin_box::testvaluechanged()
{
    Parameter testParameter(static_cast<sint>(0),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);
    Moris_Int_Spin_Box testValueChanged(nullptr, testParameter);
    testValueChanged.setObjectName("testIntSpinWorking");

    QSignalSpy setValueSpy(&testValueChanged, &Moris_Int_Spin_Box::value_changed);
    QVERIFY(setValueSpy.isValid());

    testValueChanged.setValue(15);

    QCOMPARE(setValueSpy.count(),1);
    auto spyArgs = setValueSpy.takeFirst();

    QCOMPARE(spyArgs.at(0).toString(), QString("testIntSpinWorking"));
    QCOMPARE(spyArgs.at(1).toInt(), 15);

    QCOMPARE(testParameter.get_value<sint>(),15);
}

void ut_cl_int_spin_box::testReadOnlyParameter()
{
    Parameter testParameter(static_cast<sint>(3),
                            Entry_Type::FREE,
                            "",
                            Module_Type::END_ENUM,
                            0);
    testParameter.set_value("lockedParam", static_cast<sint>(3),true);
    
    Moris_Int_Spin_Box testReadOnly(nullptr,testParameter);
    QVERIFY(testReadOnly.isReadOnly());

}