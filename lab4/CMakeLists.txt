include_directories(include)

add_executable(lab4_task1
        cli/task1.cpp
)

add_executable(lab4_task2
        cli/task2.cpp
)


#add_library(lab4_methods STATIC
#        src/ode.cpp
#        src/cli.cpp
#)
#target_include_directories(lab4_methods PUBLIC include)
#target_compile_features(lab4_methods PUBLIC cxx_std_20)
#
#foreach(task IN ITEMS 1 2)
#    add_executable(lab4_task${task} cli/task${task}.cpp)
#    target_link_libraries(lab4_task${task} PRIVATE lab4_methods)
#endforeach()
#
#if(MSVC)
#    target_compile_options(lab4_methods PRIVATE /W4 /utf-8)
#    target_compile_options(lab4_task1 PRIVATE /W4 /utf-8)
#    target_compile_options(lab4_task2 PRIVATE /W4 /utf-8)
#else()
#    target_compile_options(lab4_methods PRIVATE -Wall -Wextra -Wpedantic)
#    target_compile_options(lab4_task1 PRIVATE -Wall -Wextra -Wpedantic)
#    target_compile_options(lab4_task2 PRIVATE -Wall -Wextra -Wpedantic)
#endif()
#
#option(LAB4_BUILD_TESTS "Build numerical checks for lab4" OFF)
#if(LAB4_BUILD_TESTS)
#    enable_testing()
#    add_executable(lab4_numerical_tests tests/numerical_tests.cpp)
#    target_link_libraries(lab4_numerical_tests PRIVATE lab4_methods)
#    add_test(NAME lab4_numerical_methods COMMAND lab4_numerical_tests)
#endif()
