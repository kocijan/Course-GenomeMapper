include(FetchContent)

FetchContent_Declare(
  bioparser
  GIT_REPOSITORY https://github.com/rvaser/bioparser
  GIT_TAG 1d059386f30ffc2e2902b76cf4fb33c446a5b9bb # 3.0.12
)

FetchContent_MakeAvailable(bioparser)

add_executable(${PROJECT_NAME} crimson_mapper.cpp)
target_link_libraries(${PROJECT_NAME} bioparser)

