set(BoostVersion 1.57.0)
set(BoostMD5 1be49befbdd9a5ce9def2983ba3e7b76)

string(REGEX REPLACE "beta\\.([0-9])$" "beta\\1" BoostFolderName ${BoostVersion})
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})

ExternalProject_Add(Boost
    PREFIX Boost
    URL  http://sourceforge.net/projects/boost/files/boost/${BoostVersion}/${BoostFolderName}.tar.bz2/download
    URL_MD5 ${BoostMD5}
    CONFIGURE_COMMAND ./bootstrap.sh
                                                        --with-libraries=program_options
    BUILD_COMMAND           ./b2 install 
                                                        variant=release
                                                        link=static
                                                        cxxflags='-fPIC'
                                                        --prefix=${CMAKE_BINARY_DIR}/extern
                                                        -d 0
                                                        -j8
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

set(Boost_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
set(Boost_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/boost/ )

set(Boost_LIBRARIES -lprogram_options)
