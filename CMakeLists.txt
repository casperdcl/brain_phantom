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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")  # -Wall
  set(CMAKE_LDFLAGS "${CMAKE_LDFLAGS} -lm")
endif()

add_executable(dncat_bin ${HEADRES} ${SRCS}
  "${SRC_DIR}/dncat_output_bin.cpp")
install(TARGETS dncat_bin RUNTIME DESTINATION .)

macro(processPar PARBASE)
  add_custom_command(OUTPUT "${CMAKE_SOURCE_DIR}/${PARBASE}_1.bin"
    COMMAND dncat_bin ${PARBASE}.par ${PARBASE}
    MAIN_DEPENDENCY dncat_bin
    DEPENDS "${CMAKE_SOURCE_DIR}/${PARBASE}.par"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Running ${PARBASE}"
    VERBATIM
  )
  add_custom_target(show${PARBASE}
    COMMAND python -BO showslice.py ${PARBASE}_1.bin
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    DEPENDS "${CMAKE_SOURCE_DIR}/${PARBASE}_1.bin"
    COMMENT "Showing ${PARBASE}"
    VERBATIM
  )
endmacro()

processPar(test)
processPar(pet)
