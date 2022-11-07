# CMake generated Testfile for 
# Source directory: /home/dproksch/master_thesis/CronosCelerityBlock
# Build directory: /home/dproksch/master_thesis/CronosCelerityBlock/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tiny_run "/home/dproksch/master_thesis/CronosCelerityBlock/build/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_tiny" "SodToro5HllcRef")
set_tests_properties(tiny_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_TINY" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;201;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(tiny_compare "HDF5_DIFF_EXECUTABLE-NOTFOUND" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5")
set_tests_properties(tiny_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_TINY" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;202;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(small_run "/home/dproksch/master_thesis/CronosCelerityBlock/build/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_small" "SodToro5HllcRef")
set_tests_properties(small_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;204;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(small_compare "HDF5_DIFF_EXECUTABLE-NOTFOUND" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5")
set_tests_properties(small_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_SMALL" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;205;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(medium_run "/home/dproksch/master_thesis/CronosCelerityBlock/build/proj" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_medium" "SodToro5HllcRef")
set_tests_properties(medium_run PROPERTIES  DEPENDS "proj" FIXTURES_SETUP "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;207;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
add_test(medium_compare "HDF5_DIFF_EXECUTABLE-NOTFOUND" "--exclude-path" "/Data/version" "/home/dproksch/master_thesis/CronosCelerityBlock/build/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5" "/home/dproksch/master_thesis/CronosCelerityBlock/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5")
set_tests_properties(medium_compare PROPERTIES  FIXTURES_REQUIRED "EXEC_MEDIUM" _BACKTRACE_TRIPLES "/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;208;add_test;/home/dproksch/master_thesis/CronosCelerityBlock/CMakeLists.txt;0;")
