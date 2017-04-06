cmake_minimum_required (VERSION 3.5)
project (Bridger)

include_directories(
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}
    )


find_package( Boost 1.47 COMPONENTS system serialization thread filesystem REQUIRED)
include_directories("${Boost_INCLUDE_DIR}")

add_library(lemon STATIC IMPORTED)
set_target_properties(lemon PROPERTIES 
    IMPORTED_LOCATION /home/fin/libemon.a)


#set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
#link_directories(
#    ${PROJECT_SOURCE_DIR}/lemon
#    )

add_subdirectory(src)


#Assemble_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS) 
##Assemble_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
#Assemble_LDFLAGS =  -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
#Assemble_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)
#
#PathSearch_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS)
##PathSearch_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
#PathSearch_LDFLAGS =  -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
#PathSearch_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)
#
#get_reads_from_graph_SOURCES = utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp compatibility_graph.cpp reachability_bp_graph.cpp matching_merge.cpp get_reads_from_graph.cpp
#get_reads_from_graph_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS)
##get_reads_from_graph_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
#get_reads_from_graph_LDFLAGS = -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
#get_reads_from_graph_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)
#
#
#
