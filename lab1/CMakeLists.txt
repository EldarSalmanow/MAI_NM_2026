add_library(core
        src/vector.cpp
        src/matrix.cpp
        src/io.cpp
)

target_include_directories(core
        PRIVATE
        include
)

include_directories(include)
link_libraries(core)

add_executable(task1
        cli/task1.cpp
)
