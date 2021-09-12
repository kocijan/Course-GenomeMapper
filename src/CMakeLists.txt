include(FetchContent)

FetchContent_Declare(
  bioparser
  GIT_REPOSITORY https://github.com/rvaser/bioparser
  GIT_TAG 1d059386f30ffc2e2902b76cf4fb33c446a5b9bb # 3.0.12
)

FetchContent_MakeAvailable(bioparser)

# include_directories(../include)
# add_library(crimson_alignment_engine crimson_alignment_engine.cpp)
# add_library(crimson_minimizer_engine crimson_minimizer_engine.cpp)

add_executable(${PROJECT_NAME} crimson_mapper.cpp
${PROJECT_SOURCE_DIR}/include/crimson_alignment_engine.hpp
${PROJECT_SOURCE_DIR}/include/crimson_minimizer_engine.hpp
)

target_link_libraries(${PROJECT_NAME} PUBLIC bioparser)
target_link_libraries(${PROJECT_NAME} PUBLIC crimson_alignment_engine)
target_link_libraries(${PROJECT_NAME} PUBLIC crimson_minimizer_engine)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)