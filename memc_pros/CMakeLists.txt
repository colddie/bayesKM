cmake_minimum_required(VERSION 3.12)
 
project(meKineticRigid)
# set_target_properties(idlrtkInterface PROPERTIES LINKER_LANGUAGE CXX)

set (CMAKE_CXX_FLAGS "-fPIC")

foreach(p
    CMP0042 # CMake 3.0
    CMP0063 # CMake 3.3.2
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()


if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
    set (CMAKE_INSTALL_PREFIX "${PARENT_DIR}/bin" CACHE PATH "default install path" FORCE )
Message("CMAKE_INSTALL_PATH is self-defined!")
endif()


if( NOT CMAKE_BUILD_TYPE )
set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()

find_package(SimpleITK)
if (simpleITK_FOUND)
       include(${SimpleITK_USE_FILE})
else (simpleITK_FOUND)
	# message (FATAL_ERROR "Cannot build without SimpleITK Please set SimpleITK_DIR.")
endif (simpleITK_FOUND)

find_package (ITK)
if (ITK_FOUND)
	include ( ${ITK_USE_FILE} )
else (ITK_FOUND)
	message (FATAL_ERROR "Cannot build without ITK, Please set ITK_DIR.")
endif (ITK_FOUND)

#option(ITK_USE_GPU "GPU acceleration via OpenCL" OFF)
#mark_as_advanced(ITK_USE_GPU)

# find_package ( Elastix )
# if (ELASTIX_FOUND)
#        include ( ${ELASTIX_USE_FILE} )
# else (ELASTIX_FOUND)
# 	# message (FATAL_ERROR "Cannot build without ELASTIX, Please set ELASTIX_DIR.")
# endif (ELASTIX_FOUND)


set(ARMA_USE_FILE /home/tsun/bin/armadillo-9.400.4/install/include)
set(OPTIM_USE_FILE  /home/tsun/bin/optim-master/install/include/optim)
set(TPC_USE_PATH  /home/tsun/bin/tpcclib-master/src/v1)
set(TPC_STATIC_PATH /home/tsun/bin/tpcclib-master/build/lib/include)
set(SPLINE_USE_PATH /home/tsun/bin/spline-master/src)

option(SimpleITK_TESTING_NOSHOW "Disable the Show function in SimpleITK testing." ON)
mark_as_advanced(SimpleITK_TESTING_NOSHOW)
Message(STATUS "external libraries:" ${SimpleITK_LIBRARIES})
Message(STATUS "external libraries:" ${ITK_LIBRARIES})
Message(STATUS "external libraries:" ${ELASTIX_LIBRARIES})




# link_directories ( /storage0/home/tsun/bin/elastix-develop/install )
#link_directories ( /cm/shared/apps/hdf5/1.10.1/lib/ )
link_directories(/home/tsun/bin/tpcclib-master/build/lib)

include_directories (${TPC_USE_PATH}/libtpcmisc) 
include_directories (${TPC_USE_PATH}/libtpccurveio) 
include_directories (${TPC_USE_PATH}/libtpcmodel) 
include_directories (${TPC_USE_PATH}/libtpcsvg) 
include_directories (${TPC_USE_PATH}/libtpcimgio) 
include_directories (${TPC_USE_PATH}/libtpcimgp) 
include_directories (${TPC_USE_PATH}/libtpcmodext) 
include_directories (/home/tsun/bin/tpcclib-master/build/lib/include)

#add_executable ( meKineticRigid meKineticRigid.cxx libtpcimgp.h libtpcmisc.h tpcclibConfig.h  )  #patlak_idl.cxx
add_library(meKineticRigid SHARED meKineticRigid.cxx patlak_c.cxx meKineticRigid.h libtpcimgp.h libtpcmisc.h)
#file(GLOB ITK_LIBRARIES "/storage0/home/tsun/bin/InsightToolkit-4.13.1/bin/lib/*.so")

target_include_directories (meKineticRigid PRIVATE ${ARMA_USE_FILE} ${OPTIM_USE_FILE} ${SPLINE_USE_PATH} ${TPC_STATIC_PATH})

#Dynamic/Shared Libs
target_link_libraries(meKineticRigid ${SimpleITK_LIBRARIES} ${ITK_LIBRARIES} 
/home/tsun/bin/armadillo-9.400.4/install/lib64/libarmadillo.so 
/home/tsun/bin/optim-master/install/lib/liboptim.so
gomp)

target_link_libraries(meKineticRigid 
/home/tsun/bin/tpcclib-master/build/lib/libtpcimgp.a
/home/tsun/bin/tpcclib-master/build/lib/libtpcimgio.a
/home/tsun/bin/tpcclib-master/build/lib/libtpccurveio.a
/home/tsun/bin/tpcclib-master/build/lib/libtpcmisc.a
/home/tsun/bin/tpcclib-master/build/lib/libtpcmodel.a
/home/tsun/bin/tpcclib-master/build/lib/libtpcmodext.a
/home/tsun/bin/tpcclib-master/build/lib/libtpcsvg.a)
#tpcimgp_shared tpcimgio_shared tpccurveio_shared tpcmodext_shared
#tpcmodel_shared tpccurveio_shared tpcmisc_shared)




#install(TARGETS rawTomha 
#         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
#install(TARGETS  rawTomhalib 
#         ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX})		 
install(TARGETS  meKineticRigid
         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}         
         LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX})     ## was ARCHIVE
		 
