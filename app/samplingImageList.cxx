/*=========================================================================
 Program:   ORFEO Toolbox
 Language:  C++
 Date:      $Date$
 Version:   $Revision$


 Copyright (c) Centre National d'Etudes Spatiales. All rights reserved.
 See OTBCopyright.txt for details.


 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "otbImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbOGRDataSourceWrapper.h"
#include "otbMultiChannelExtractROI.h"
#include "otbExtractROI.h"
#include "otbOGRFeatureWrapper.h"
#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"
#include "otbWrapperOutputFilenameParameter.h"
#include "itkSmartPointer.h"
#include "itkPreOrderTreeIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkListSample.h"
#include "vcl_algorithm.h"
#include "otbMultiToMonoChannelExtractROI.h"
#include <sstream>
#include <iterator>   

namespace otb
{

namespace Wrapper
{
  
template <typename T>
std::string to_string(T value)
{
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}   
  
class SamplingImageList : public Application
{
public:
  typedef SamplingImageList Self;
  typedef itk::SmartPointer<Self> Pointer; 
  
  typedef FloatVectorImageType                 ImageType;  
  typedef ImageType::InternalPixelType         ImagePixelType;
  typedef UInt32ImageType                      LabelImageType;
  typedef LabelImageType::InternalPixelType    LabelImagePixelType; 
  typedef itk::ImageRegionIterator<ImageType>  IteratorType;
  ImageType::IndexType IndexType;
  
  typename ImageType::RegionType            polygonRegion;
  
  typedef otb::MultiChannelExtractROI<ImagePixelType, ImagePixelType>  ExtractROIFilterType;
  
  itkNewMacro(Self);

  itkTypeMacro(SamplingImageList, otb::Application);
  
private:
  void DoInit()
  {
    SetName("SamplingImageList");
    SetDescription("This application generates a .txt file from an image and a vectoriel file. ");
    
    AddParameter(ParameterType_InputImageList, "in", "Input Image");    
    
    AddParameter(ParameterType_InputFilename, "shp", "Vectoriel File");    
    AddParameter(ParameterType_OutputFilename, "out", "Output Text");    
    AddParameter(ParameterType_OutputFilename, "v", "Verification Mask");    
    AddParameter(ParameterType_String, "cfield", "Field of class");
    
    AddParameter(ParameterType_Choice, "mode", "Mode of sampling");
    AddChoice("mode.exhaustive", "Exhaustive Sampling");
    AddChoice("mode.random", "Random Sampling");
    AddChoice("mode.randomequally", "Random sampling equally distributed according to the classes size");
    AddChoice("mode.periodic", "Periodic sampling, in all the polygons");
    AddChoice("mode.periodicrandom", "Periodic sampling, in all the polygons, randomly shifted");
        
    AddParameter(ParameterType_Int, "samples", "Number of samples per classes");
    SetDefaultParameterInt("samples", 2000);
    MandatoryOff("samples"); 
    
    AddParameter(ParameterType_Int, "tiles", "Size of square tiles");
    SetDefaultParameterInt("tiles", 500);
    MandatoryOff("tiles");
    
    AddParameter(ParameterType_Int, "nd", "NoData value");
    SetDefaultParameterInt("nd", 0);
    MandatoryOff("nd"); 
    
    AddParameter(ParameterType_Int, "rand", "Seed value for Mersenne Twister Random Generator");
    MandatoryOff("rand"); 
     
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  {  
    FloatVectorImageListType::Pointer inList = this->GetParameterImageList("in");

    if( inList->Size() == 0 )
    {
      itkExceptionMacro("No input Image set...");
    }
    
    for( unsigned int i=0; i<inList->Size(); i++ )
    {
      inList->GetNthElement(i)->UpdateOutputInformation();
    }
    
    //Output text file
    std::ofstream myfile;
    myfile.open (GetParameterString("out").c_str());
    
    //Input shape file
    otb::ogr::DataSource::Pointer vectorData = otb::ogr::DataSource::New(GetParameterString("shp").c_str(), otb::ogr::DataSource::Modes::Read);
    
    //Output shape file
    std::vector<std::string> options;
    std::string projRef = inList->GetNthElement(0)->GetProjectionRef();
    OGRSpatialReference oSRS(projRef.c_str());
    
    otb::ogr::DataSource::Pointer ogrDS;
    ogrDS = otb::ogr::DataSource::New(GetParameterString("v").c_str(), otb::ogr::DataSource::Modes::Overwrite);
    otb::ogr::Layer layer(NULL, false);
    std::string layername = itksys::SystemTools::GetFilenameName(GetParameterString("v").c_str());
    std::string extension = itksys::SystemTools::GetFilenameLastExtension(GetParameterString("v").c_str());
    layername = layername.substr(0,layername.size()-(extension.size()));
    layer = ogrDS->CreateLayer(layername, &oSRS, wkbPoint, options);
    
    //No data value for pixels value
    int noDataValue = GetParameterInt("nd");
        
    //Dimension of a side of the square tiles
    int sizeTiles = GetParameterInt("tiles");
    
    //Seed value for Mersenne Twister Random Generator 
    int seed = GetParameterInt("rand");
    
    //Mode recuperation from the parameter
    const std::string samplingMode = GetParameterString("mode");
    
    //Number of elements in each pixels
    unsigned int nbComponents = inList->GetNthElement(0)->GetNumberOfComponentsPerPixel();
    
    //Initialisation of the output shape file
    otb::ogr::Layer preFiltered = vectorData->GetLayer(0);    
    otb::ogr::Feature preFeature = preFiltered.GetFeature(0);
        
    for(unsigned int comp = 0; comp<nbComponents; ++comp)
    {  
      std::ostringstream fieldoss;
      fieldoss<<"b"<<comp;
      OGRFieldDefn field(fieldoss.str().c_str(), OFTReal);
      layer.CreateField(field, true);
    }  
    
    for(int comp = 0; comp < preFeature.ogr().GetFieldCount(); ++comp)
    {  
      std::ostringstream fieldoss;
      fieldoss<<comp + nbComponents;
      OGRFieldDefn field(preFeature.ogr().GetFieldDefnRef(comp)->GetNameRef(), preFeature.ogr().GetFieldDefnRef(comp)->GetType());
      layer.CreateField(field, true);  
    }
    
    unsigned long sizeTilesX = sizeTiles;
    unsigned long sizeTilesY = sizeTiles;
    std::map<unsigned int, unsigned long> sizeImageX;
    std::map<unsigned int, unsigned long> sizeImageY;    
    std::map<unsigned int, unsigned long> nbTilesX;
    std::map<unsigned int, unsigned long> nbTilesY;
    
    for( unsigned int i=0; i<inList->Size(); i++ )
    {
      sizeImageX[i] = inList->GetNthElement(i)->GetLargestPossibleRegion().GetSize()[0];
      sizeImageY[i] = inList->GetNthElement(i)->GetLargestPossibleRegion().GetSize()[1];
      nbTilesX[i] = sizeImageX[i]/sizeTilesX + (sizeImageX[i]%sizeTilesX > 0 ? 1 : 0);
      nbTilesY[i] = sizeImageY[i]/sizeTilesY + (sizeImageY[i]%sizeTilesY > 0 ? 1 : 0);
    }
    
    //Number of pixels in all the polygons
    int nbPixelsGlobal = 0; 
    
    //Number of pixels in each classes
    std::map<int, int> elmtsInClass;
    
    //Number of pixels in each polygons
    std::map<unsigned long, int> polygon;
    //Counter of pixels in the current polygon
    std::map<unsigned long, int> counterPixelsInPolygon;
    //Shifted counter of pixels in the current polygon
    std::map<unsigned long, int> counterPixelsInPolygonShifted;
    //RandomPosition of a sampled pixel for each polygons
    std::map<unsigned long, int> randomPositionInPolygon;

    //Varibles to build the progression bar
    int stepsProgression = 0;
    int currentProgression;
    
    otbAppLogINFO(<< "Computing the number of pixels for each polygons and classes" << std::endl);
    
    //Initialisation of the random generator
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize();
    generator->SetSeed(seed);
    
    int polyForced = 0;
    
    
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::SamplingImageList)