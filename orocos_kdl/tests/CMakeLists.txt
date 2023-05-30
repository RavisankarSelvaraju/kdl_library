 MESSAGE(STATUS "test enabled")
 set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
 INCLUDE_DIRECTORIES(${PROJ_SOURCE_DIR}/src ${PROJ_SOURCE_DIR}/models ${PROJ_BINARY_DIR}/src)

#  add_executable(solvertest solvertest.cpp )
 add_executable(test_vereshchagin test_vereshchagin.cpp )
 TARGET_LINK_LIBRARIES(test_vereshchagin PRIVATE orocos-kdl cppunit)
