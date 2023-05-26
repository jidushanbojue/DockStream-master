#cmake_minimum_required(VERSION 3.18)
#project(AutoDock-Vina-1.2.3)
##set(CMAKE_CXX_STANDARD 14)
##include_directories(src/lib)
#file(GLOB_RECURSE VINA_LIB_SRC src/lib/*.cpp)
#add_library(vina_lib STATIC ${VINA_LIB_SRC})
#
#add_executable(vina src/main/main.cpp)
##add_executable(vina_split src/split/main.cpp)
#target_link_libraries(vina vina_lib)

cmake_minimum_required(VERSION 3.18)
project(AutoDock-Vina-1.2.3)
#set(CMAKE_CXX_STANDARD 14)
include_directories(src/lib)
#add_executable(vina src/main/main.cpp src/lib/ad4cache.cpp src/lib/conf_independent.cpp  src/lib/grid.cpp
#        src/lib/monte_carlo.cpp src/lib/non_cache.cpp src/lib/parallel_progress.cpp src/lib/quasi_newton.cpp
#        src/lib/random.cpp src/lib/utils.cpp src/lib/cache.cpp src/lib/coords.cpp src/lib/model.cpp
#        src/lib/mutate.cpp src/lib/parallel_mc.cpp src/lib/parse_pdbqt.cpp src/lib/quaternion.cpp
#        src/lib/szv_grid.cpp src/lib/vina.cpp)

add_executable(vina src/main/main.cpp src/lib/utils.cpp src/lib/vina.cpp)