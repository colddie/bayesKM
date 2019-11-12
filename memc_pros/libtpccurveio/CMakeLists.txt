project(libtpccurveio)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpccurveio SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpccurveio PRIVATE ../include)