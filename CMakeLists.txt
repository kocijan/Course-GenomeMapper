cmake_minimum_required(VERSION 3.21.0)
project(crimson_mapper VERSION 0.1.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

configure_file(include/${PROJECT_NAME}Config.h.in include/${PROJECT_NAME}Config.h)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(include)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -O0 -g")
  set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -O0 -g")
endif()

# include(CTest)
enable_testing()
# add_test(NAME vtest COMMAND ${PROJECT_NAME} --version)
# add_test(NAME ecoli1 COMMAND ${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/test/data/GCF_000005845.2_ASM584v2_genomic.fna ${PROJECT_SOURCE_DIR}/test/data/MAP006-1_2D_pass.fasta)
# add_test(NAME ecoli12 COMMAND ${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/test/data/GCF_000005845.2_ASM584v2_genomic.fna ${PROJECT_SOURCE_DIR}/test/data/MAP006-1_2D_pass.fasta ${PROJECT_SOURCE_DIR}/test/data/MAP006-2_2D_pass.fasta)

target_compile_options(${PROJECT_NAME} PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wconversion>
)
target_compile_options(crimson_alignment_engine PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wconversion>
)
