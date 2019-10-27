project(jsoncpp) # the name of your project

cmake_minimum_required(VERSION 2.4.0)

# tell cmake to process CMakeLists.txt in that subdirectory
add_subdirectory(src/lib_json)
add_subdirectory(include/json)


