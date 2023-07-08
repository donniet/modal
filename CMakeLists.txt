cmake_minimum_required(VERSION 3.10)

project(modal)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(TBB REQUIRED)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

enable_testing()

add_executable(
  ${CMAKE_PROJECT_NAME}_test
  src/main_test.cpp
)
target_link_libraries(
  ${CMAKE_PROJECT_NAME}_test PRIVATE
  gtest_main
  TBB::tbb
)

include(GoogleTest)
gtest_discover_tests(${CMAKE_PROJECT_NAME}_test)



add_executable(${CMAKE_PROJECT_NAME} src/main.cpp) 
