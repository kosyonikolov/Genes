add_executable(CameraTest src/cameraTest.cpp)

find_package(OpenCV REQUIRED)
target_link_libraries(CameraTest PRIVATE ${OpenCV_LIBS})
target_link_libraries(CameraTest PRIVATE Genetic)