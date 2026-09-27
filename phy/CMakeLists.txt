project(phy)

file(GLOB_RECURSE SOURCES *.cpp)

add_executable(phy ${SOURCES})

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBIIO REQUIRED IMPORTED_TARGET libiio)

target_compile_features(phy PRIVATE cxx_std_17)

target_link_libraries(phy
    PRIVATE
        PkgConfig::LIBIIO
)
