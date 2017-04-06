add_executable(Assemble utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp assemble.cpp)
target_link_libraries(Assemble ${Boost_LIBRARIES})

add_executable(PathSearch utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp compatibility_graph.cpp reachability_bp_graph.cpp matching_merge.cpp path_search.cpp)
target_link_libraries(PathSearch ${Boost_LIBRARIES} lemon)
    #${PROJECT_SOURCE_DIR}/libemon.a)

add_executable(get_reads_from_graph utility.cpp kmerhash.cpp common.cpp loadreads.cpp splicing_graph.cpp compatibility_graph.cpp reachability_bp_graph.cpp matching_merge.cpp get_reads_from_graph.cpp)
target_link_libraries(get_reads_from_graph ${Boost_LIBRARIES})
