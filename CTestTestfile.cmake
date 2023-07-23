# CMake generated Testfile for 
# Source directory: /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity
# Build directory: /home/dproksch/Documents/master/master_thesis/repos/CronosCelerity
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.

# add_test(tiny_run "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_tiny" "SodToro5HllcRef")
# set_tests_properties(tiny_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_TINY" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;293;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(tiny_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5")
# set_tests_properties(tiny_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_TINY" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;294;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(small_run "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_small" "SodToro5HllcRef")
# set_tests_properties(small_run PROPERTIES  DEPENDS "proj" FAIL_REGULAR_EXPRESSION "Error;error;ERROR" FIXTURES_SETUP "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;296;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(small_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5")
# set_tests_properties(small_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;297;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(medium_run "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_medium" "SodToro5HllcRef")
# set_tests_properties(medium_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;299;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(medium_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5")
# set_tests_properties(medium_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;300;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(new_test_run "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_new_test" "SodToro5HllcRef")
# set_tests_properties(new_test_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_NEW_TEST" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;299;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")

# add_test(tiny_run "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_tiny" "SedovCartesian_hllc")
# set_tests_properties(tiny_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_TINY" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;235;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
# add_test(tiny_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step1.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step1.h5")
# set_tests_properties(tiny_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;297;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(small_run "mpirun -n 2" "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_small" "SedovCartesian_hllc")
set_tests_properties(small_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;236;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(small_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5")
set_tests_properties(small_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;298;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(medium_run "mpirun -n 2" "/home/dproksch/master_thesis/CronosCelerityBlock/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_medium" "SedovCartesian_hllc")
set_tests_properties(medium_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;237;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(medium_compare "/usr/bin/h5diff" "-v" "-d 1.5e-11" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step29.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step29.h5")
set_tests_properties(medium_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;298;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")

subdirs("external/CronosNumLib")
