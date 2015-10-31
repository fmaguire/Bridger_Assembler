AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I ../m4
# Libtool dependencies
LIBTOOL_DEPS = @LIBTOOL_DEPS@

bin_PROGRAMS = Assemble PathSearch get_reads_from_graph

Assemble_SOURCES = utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp assemble.cpp utility.h kmerhash.h common.h loadreads.h splicing_graph.h assemble.h
Assemble_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS) 
#Assemble_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
Assemble_LDFLAGS =  -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
Assemble_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)

PathSearch_SOURCES = utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp compatibility_graph.cpp reachability_bp_graph.cpp matching_merge.cpp path_search.cpp
PathSearch_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS)
#PathSearch_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
PathSearch_LDFLAGS =  -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
PathSearch_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)

get_reads_from_graph_SOURCES = utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp compatibility_graph.cpp reachability_bp_graph.cpp matching_merge.cpp get_reads_from_graph.cpp
get_reads_from_graph_CPPFLAGS =  -g -Wno-deprecated -static $(CPPFLAGS)
#get_reads_from_graph_LDFLAGS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_LDFLAGS)
get_reads_from_graph_LDFLAGS = -lboost_serialization -lboost_system -lboost_filesystem $(BOOST_LDFLAGS)
get_reads_from_graph_LIBS = $(BOOST_THREAD_LIBS) $(BOOST_SERIALIZATION_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS)
