CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(oofdtd)

SET(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
#SET(CMAKE_CXX_COMPILER /opt/apps/intel15/mvapich2/2.1/bin/mpicxx)

SET(CMAKE_CXX_FLAGS "-std=c++11 -O2 -w -Wall -Werror")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE(STATUS "  *************** ")
FOREACH(item ${CMAKE_MODULE_PATH})
    MESSAGE(STATUS "  path: " ${item})
ENDFOREACH()
MESSAGE(STATUS "  *************** ")

# package config
FIND_PACKAGE( PkgConfig REQUIRED )

FILE(GLOB oofdtd_SRC "src/*.cpp")

SET(INCL_FOLDER "include")
INCLUDE_DIRECTORIES(${INCL_FOLDER})

ADD_EXECUTABLE(oofdtd main.cpp ${oofdtd_SRC})

FOREACH(item ${})
    MESSAGE(STATUS "  including: " ${item})
ENDFOREACH()
FOREACH(item ${ALL_LIBRARIES})
    MESSAGE(STATUS "  using lib: " ${item})
ENDFOREACH()

TARGET_LINK_LIBRARIES( oofdtd ${ALL_LIBRARIES})
