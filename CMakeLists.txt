cmake_minimum_required(VERSION 3.1.0)

set(tid 06T)

project(ceng391_${tid} CXX)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SRCS image.cc util.cc homography.cc)
set(HDRS image.hpp
  types.hpp
  keypoint.hpp
  descriptor.hpp
  brief_data.hpp
  util.hpp
  homography.hpp)

find_package(BLAS)
find_package(LAPACK)

add_executable(image-test ${SRCS} image_test.cc)
target_link_libraries(image-test ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})

add_custom_target(archive rm ceng391_${tid}/*~ -vf && tar --exclude "ceng391_${tid}/build" --exclude "compile_commands.json" --exclude ".clangd" -czvf ceng391_${tid}.tar.gz ceng391_${tid}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
  COMMENT "Pack and compress files for distribution")

add_custom_target(edit emacsclient -n ${SRCS} ${HDRS} image_test.cc CMakeLists.txt
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Open source and header files using Emacs server")

add_custom_target(test-valgrind valgrind --leak-check=full --show-leak-kinds=all ./image-test
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build
  DEPENDS image-test
  COMMENT "Check for memory leaks using Valgrind")
