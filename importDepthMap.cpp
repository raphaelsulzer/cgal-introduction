#include <cgal_typedefs.h>
#include "tiffio.h"


void readTiff(std::string path)
{
    TIFF* tif = TIFFOpen(path.c_str(), "r");

//    if(tif){
//        uint32 w, h;
//        size_t npixels;
//        double* raster;

//        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
//        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
//        npixels = w * h;
//        raster = (double*) _TIFFmalloc(npixels * sizeof (uint32));
//        TIFFReadRGBAImage(tif, w, h, raster, 0);
//        auto data1 = raster[0,0];
////        if (raster != NULL) {
////           if (TIFFReadRGBAImage(tif, w, h, raster, 0)) {
////           ...process raster data...
////           }
////           _TIFFfree(raster);
////        }
//        TIFFClose(tif);
//    }


}
