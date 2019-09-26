
# Make sure the compiler can find include files from our libraries. 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcmisc) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpccurveio) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcmodel) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcsvg) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcimgio) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcimgp) 
include_directories (${CMAKE_SOURCE_DIR}/v1/libtpcmodext) 

# Make sure the linker can find the libraries once they are built. 
link_directories (${CMAKE_BINARY_DIR}/v1/libtpcmisc) 
link_directories (${CMAKE_BINARY_DIR}/v1/libtpccurveio) 
link_directories (${CMAKE_BINARY_DIR}/v1/libtpcmodel) 
link_directories (${CMAKE_BINARY_DIR}/v1/libtpcsvg) 
link_directories (${CMAKE_BINARY_DIR}/v1/libtpcmodext) 

# Add executable that is built from the source file(s) 
# The extensions are automatically found. 
add_executable (logan logan.c) 
add_executable (patlak patlak.c) 
add_executable (regfur regfur.c) 
add_executable (imgfur imgfur.c) 
add_executable (imgki imgki.c) 
add_executable (imgdv imgdv.c) 
# add_executable (patlak_idl patlak_idl.c) 

add_library (mtga_idl SHARED 
    patlak_idl.c logan_idl.c regfur_idl.c mrtm_idl.c simPatlak.c simLogan.c simPatlak_idl.c simLogan_idl.c 
    tcm2_idl.c tcm2_reverse_idl.c srtm_idl.c sim2cm.c
)
set_property(TARGET mtga_idl PROPERTY POSITION_INDEPENDENT_CODE ON)

# Link the executable to the libraries. 
target_link_libraries (logan  tpcmodext tpcsvg tpcmodel tpccurveio tpcmisc m)
target_link_libraries (patlak tpcmodext tpcmodel tpccurveio tpcmisc m)
target_link_libraries (regfur tpcmodext tpcmodel tpccurveio tpcmisc m)
target_link_libraries (imgfur tpcmodext tpcimgp tpcimgio tpcmodel tpccurveio tpcmisc m)
target_link_libraries (imgki tpcmodext tpcimgp tpcimgio tpcmodel tpccurveio tpcmisc m)
target_link_libraries (imgdv tpcmodext tpcimgp tpcimgio tpcmodel tpccurveio tpcmisc m)
target_link_libraries ( 
  mtga_idl
  tpcmodext tpcmodel tpccurveio tpcmisc m
)


# Install the executable(s)
install(
  PROGRAMS
  ${CMAKE_CURRENT_BINARY_DIR}/logan${CMAKE_EXECUTABLE_SUFFIX}
  ${CMAKE_CURRENT_BINARY_DIR}/patlak${CMAKE_EXECUTABLE_SUFFIX}
  ${CMAKE_CURRENT_BINARY_DIR}/regfur${CMAKE_EXECUTABLE_SUFFIX}
  ${CMAKE_CURRENT_BINARY_DIR}/imgfur${CMAKE_EXECUTABLE_SUFFIX}
  ${CMAKE_CURRENT_BINARY_DIR}/imgki${CMAKE_EXECUTABLE_SUFFIX}
  ${CMAKE_CURRENT_BINARY_DIR}/imgdv${CMAKE_EXECUTABLE_SUFFIX}
  DESTINATION bin
  COMPONENT applications
)

install (
  TARGETS mtga_idl
  RUNTIME
  LIBRARY
  #ARCHIVE 
  DESTINATION bin
  COMPONENT libraries
)


# Copy test data folder
# Permissions are set to avoid 'Permission denied' error on OSX
file(
  COPY test 
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
  PATTERN "*.sh" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
    GROUP_READ GROUP_WRITE GROUP_EXECUTE
    WORLD_READ WORLD_EXECUTE
)

#
# Simple tests
#

# Does does the usage message work

add_test (loganUsage logan "--help")
set_tests_properties (loganUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)
add_test (patlakUsage patlak "--help")
set_tests_properties (patlakUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)
add_test (regfurUsage regfur "--help")
set_tests_properties (regfurUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)
add_test (imgfurUsage imgfur "--help")
set_tests_properties (imgfurUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)
add_test (imgkiUsage imgki "--help")
set_tests_properties (imgkiUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)
add_test (imgdvUsage imgdv "--help")
set_tests_properties (imgdvUsage
  PROPERTIES PASS_REGULAR_EXPRESSION "Usage: *"
)

#
# Run test scripts
# 
add_test (
  NAME loganTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/logan
  COMMAND bash "./test_logan.sh"
)
add_test (
  NAME patlakTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/patlak
  COMMAND bash "./test_patlak.sh"
)
add_test (
  NAME regfurTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/regfur
  COMMAND bash "./test_regfur.sh"
)
add_test (
  NAME imgfurTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/imgfur
  COMMAND bash "./test_imgfur.sh"
)
add_test (
  NAME imgkiTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/imgki
  COMMAND bash "./test_imgki.sh"
)
add_test (
  NAME imgdvTests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test/imgdv
  COMMAND bash "./test_imgdv.sh"
)

