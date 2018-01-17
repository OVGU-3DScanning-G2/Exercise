project(3dS_01_1)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_AUTOMOC ON)

# inlude OpenMP if installed:
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(INCLUDE
    include/SVD.h
    include/Algorithms.h
    include/KDTree.h
    include/GLcamera.h
    include/Point3d.h
    include/Shader.h)

set(SOURCE
    source/SVD.cpp
    source/Algorithms.cpp
    source/KDTree.cpp
    source/GLcamera.cpp
    source/OpenGL_Application.cpp
    source/Point3d.cpp
    source/Shader.cpp)

aux_source_directory(GLFW EXTERN)
set (CMAKE_CXX_STANDARD 11)
add_executable(${PROJECT_NAME}
                ${SOURCE}
#                ${EXTERN}
                ${INCLUDE}
)


#find_package(PkgConfig REQUIRED)
#pkg_search_module(GLFW REQUIRED glfw3)

find_package(glfw3 REQUIRED)
#find_package(Qt5Widgets REQUIRED) # Core and Gui are also included
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)



#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#include_directories( ${OPENGL_INCLUDE_DIRS}
##                     ${GLUT_INCLUDE_DIRS}
#                     ${Qt5Widgets_INCLUDE_DIRS}
##                     ${CMAKE_CURRENT_SOURCE_DIR}
#                     ${Qt5Core_INCLUDE_DIRS}
#                    )


target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY}  glfw)
