file(COPY ${CMAKE_CURRENT_LIST_DIR}/reference DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(_ymw16_tests ymw16_tests.cc)

foreach (_t ${_ymw16_tests})
  add_executable(${_t}_exe ${_t} ${GTEST_LIB_SOURCES} ${GTEST_MAIN_SOURCES})
  target_include_directories(${_t}_exe PRIVATE ${ALL_INCLUDE_DIR} ${GTEST_INCLUDE_DIRS})
  target_link_libraries(${_t}_exe ${CMAKE_THREAD_LIBS_INIT} ${ALL_LIBRARIES} hammurabi)
  add_test(${_t} ${_t}_exe)
endforeach()