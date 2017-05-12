project(brain_phantom)
cmake_minimum_required(VERSION 2.8)

include_directories("${CMAKE_SOURCE_DIR}/include")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
file(GLOB ALL_SRCS "${SRC_DIR}/*.cpp")
set(SRCS ${ALL_SRCS})
list(REMOVE_ITEM SRCS "${SRC_DIR}/dncat_output_im.cpp")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} /Ox /Oi /Ot /GT /GL /MP")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O6 -g -Wall -std=c++11")
  set(CMAKE_LDFLAGS "${CMAKE_LDFLAGS} -lm")
endif()

add_executable(dncat_bin ${HEADRES} ${SRCS}
  "${SRC_DIR}/dncat_output_bin.cpp")
install(TARGETS dncat_bin RUNTIME DESTINATION .)

add_custom_target(run
  COMMAND dncat_bin test.par test
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  DEPENDS dncat_bin test.par
  COMMENT "Running test"
  VERBATIM
)
# python -BO showslice.py test_1.bin
