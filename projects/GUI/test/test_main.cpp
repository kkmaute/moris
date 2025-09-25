#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include "ut_gui_startup.hpp"
#include "ut_cl_line_edit.hpp"
#include "ut_cl_int_spin_box.hpp"
#include "ut_cl_double_spin_box.hpp"
#include "ut_cl_combo_box.hpp"
#include "ut_cl_bool_combo_box.hpp"

#include <QApplication>
#include <QTest>

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"



moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

TEST_CASE("GUI startup launches correctly", "[QTest]"){
    ut_gui_startup  testGS; 

    int qtResult = QTest::qExec(&testGS, 0, nullptr);
    REQUIRE(qtResult == 0);
}
TEST_CASE("Line-edit widget behavior is correct", "[QTest]"){
    ut_cl_line_edit testLE; 
 
    int qtResult = QTest::qExec(&testLE, 0, nullptr);
    REQUIRE(qtResult == 0);
}
TEST_CASE("Spin-box widget behavior is correct", "[QTest]"){
    ut_cl_int_spin_box testISB;

    int qtResult = QTest::qExec(&testISB, 0, nullptr);
    REQUIRE(qtResult == 0);
}
TEST_CASE("Double spin-box widget behavior is correct", "[QTest]"){
    ut_cl_double_spin_box testDSB;

    int qtResult = QTest::qExec(&testDSB, 0, nullptr);
    REQUIRE(qtResult == 0);
}
TEST_CASE("Combo-box widget behavior is correct", "[QTest]"){
    ut_cl_combo_box testCB;

    int qtResult = QTest::qExec(&testCB, 0, nullptr);
    REQUIRE(qtResult == 0);
}
TEST_CASE("Bool combo-box widget behavior is correct", "[QTest]"){
    ut_cl_bool_combo_box testBCB;

    int qtResult = QTest::qExec(&testBCB, 0, nullptr);
    REQUIRE(qtResult == 0);
}

int
main(
    int argc,
    char * argv[] )
{
    std::cout << "GUI test made" << std::endl;

    gMorisComm.initialize(&argc, &argv);

    gLogger.initialize(2);

    QApplication app(argc, argv);

    int result = Catch::Session().run(argc, argv);
    
    // // put your tests in here:
    // ut_gui_startup  testGS; 
    // ut_cl_line_edit testLE;

    // result |= QTest::qExec(&testGS, argc, argv);
    // result |= QTest::qExec(&testLE, argc, argv);

    // result |= QTest::qExec(new ut_file_dialogs, argc, argv);

    gMorisComm.finalize();
    return result;
}

