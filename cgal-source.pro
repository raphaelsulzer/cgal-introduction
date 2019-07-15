TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    fileIO.cpp \
    rayTracing.cpp \
    surfaceRecon.cpp \
    optimization.cpp \
    meshProcessing.cpp \
    pointSetProcessing.cpp \
    plyDefinition.cpp \
    colmapPLY.cpp \
    rPLY.c \
    main.cpp \
    meshPLY.cpp

QMAKE_CXXFLAGS += -std=c++11

# surpress weird warnings
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-variable -Wno-unused-parameter

### UBUNTU
unix:!macx{
    INCLUDEPATH +=  /usr/local/include \
                    /usr/include/eigen3 \
                    /usr/lib/gco-v3.0-master \
                    /home/raphael/PhD_local/cpp/rPLY

    LIBS += \
    -L/usr/local/lib/ -lCGAL \
    -L/usr/lib/ -lgmp \
    -L/usr/lib/ -lmpfr \
    -L/usr/lib/gco-v3.0-master/build/ -lgco
#    -L/usr/lib/librply-master/bin/ -lrply
}

### MAC OS
macx{
    INCLUDEPATH +=  /usr/local/Cellar/cgal/4.13/include \ # CGAL
                    /usr/local/Cellar/gmp/6.1.2_2/include \ # GMP
                    /usr/local/Cellar/mpfr/4.0.1/include \ # MPFR
                    /usr/local/Cellar/boost/1.68.0_1/include \ # BOOST
                    /usr/local/gco-v3.0-master \ #GCoptimization
                    /usr/local/eigen3 \ #GCoptimization
                    /usr/local/Cellar/pcl/1.8.0_7/include/pcl-1.8 \ # pcl
                    /usr/local/Cellar/flann/1.9.1_4/include \ # also for pcl
                    /usr/local/ply_file_reader

    LIBS += -L/usr/local/lib/ -lgmp \
            -L/usr/local/lib/ -lmpfr \
            -L/usr/local/lib/ -lCGAL \
            -L/usr/local/gco-v3.0-master/build/ -lgco \
}




HEADERS += \
    cgal_typedefs.h
