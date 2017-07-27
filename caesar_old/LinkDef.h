
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedef;

#pragma link C++ class ImgFITSReader+;
#pragma link C++ class vector<ImgFITSReader>+;
#pragma link C++ class vector<ImgFITSReader*>+;

#pragma link C++ class FITSReader+;
#pragma link C++ class vector<FITSReader>+;
#pragma link C++ class vector<FITSReader*>+;

#pragma link C++ class Img+;
#pragma link C++ class vector<Img>+;
#pragma link C++ class vector<Img*>+;

#pragma link C++ struct Img::MetaData+;
#pragma link C++ struct Img::StatsData+;
#pragma link C++ class std::vector<Img::StatsData>+;

#pragma link C++ struct Img::BkgData+;
#pragma link C++ struct Img::BkgData*+;
#pragma link C++ class std::vector<Img::BkgData>+;
#pragma link C++ class std::vector<Img::BkgData*>+;


#pragma link C++ class BkgFinder+;
#pragma link C++ class EllipseUtils+;
#pragma link C++ class ZernikeMoments+;

#pragma link C++ class Source+;
#pragma link C++ class vector<Source>+;
#pragma link C++ class vector<Source*>+;

#pragma link C++ class Contour+;
#pragma link C++ class vector<Contour>+;
#pragma link C++ class vector<Contour*>+;
#pragma link C++ class cv::Point2f+;
#pragma link C++ class cv::Point+;
#pragma link C++ class vector<cv::Point>+;
#pragma link C++ class cv::Moments;

#pragma link C++ class Region+;
#pragma link C++ class vector<Region>+;
#pragma link C++ class vector<Region*>+;

#pragma link C++ class ChanVeseSegmentation+;

#pragma link C++ class SourceFitter+;
#pragma link C++ class SLICUtils+;
#pragma link C++ class SLICSegmenter+;
#pragma link C++ class HClust+;
#pragma link C++ class OutlierDetector+;
#pragma link C++ class Interpolator+;
#pragma link C++ class BkgFinder+;


/*
#pragma link C++ class WorldCoor+;
#pragma link C++ class poly+;
#pragma link C++ class wcsprm+;
#pragma link C++ class linprm+;
#pragma link C++ class celprm+;
#pragma link C++ class prjprm+;
#pragma link C++ class IRAFsurface+;
#pragma link C++ class Distort+;
*/
#endif

