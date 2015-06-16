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
#include "otbStatisticsXMLFileWriter.h"
#include "otbStatisticsXMLFileReader.h"
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
    SetDescription("This application sample pixels from an image and shape file with instruction from .xml file.");
    
    AddParameter(ParameterType_InputImage, "in", "Input Image List");    
    AddParameter(ParameterType_InputFilename, "xml", "XML Analysis File");
    AddParameter(ParameterType_InputFilename, "shp", "Vectoriel File"); 
    AddParameter(ParameterType_OutputFilename, "v", "Verification Mask");
    AddParameter(ParameterType_OutputFilename, "out", "Output Text");    
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
    SetDefaultParameterInt("tiles", 200);
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
    ImageType::Pointer image = GetParameterImage("in");
    image->UpdateOutputInformation(); 
    
    //Output text file
    std::ofstream myfile;
    myfile.open (GetParameterString("out").c_str());
    
    //Input shape file
    otb::ogr::DataSource::Pointer vectorData = otb::ogr::DataSource::New(GetParameterString("shp").c_str(), otb::ogr::DataSource::Modes::Read);
    
    //Output shape file
    std::vector<std::string> options;
    std::string projRef = image->GetProjectionRef();
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
    unsigned int nbComponents = image->GetNumberOfComponentsPerPixel();
    
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
        
    //Initialisation of the random generator
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize();
    generator->SetSeed(seed);
    
    
    
    int polyForced = 0;
              
    otbAppLogINFO(<< "Sampling pixels with the sampling mode : " << samplingMode << std::endl);
    
    //Tiling
    unsigned long sizeTilesX = sizeTiles;
    unsigned long sizeTilesY = sizeTiles;
    unsigned long sizeImageX = image->GetLargestPossibleRegion().GetSize()[0];
    unsigned long sizeImageY = image->GetLargestPossibleRegion().GetSize()[1];    
    unsigned int nbTilesX    = sizeImageX/sizeTilesX + (sizeImageX%sizeTilesX > 0 ? 1 : 0);
    unsigned int nbTilesY    = sizeImageY/sizeTilesY + (sizeImageY%sizeTilesY > 0 ? 1 : 0);   
     
    //Progression bar re-initialisation
    stepsProgression = 0;    
      
    TiXmlDocument doc(GetParameterString("xml").c_str());
    if(!doc.LoadFile())
    {
      std::cout << "le DOC n'existe pas" << std::endl;
    }
    TiXmlElement *elem = doc.FirstChildElement()->FirstChildElement();
    if(!elem)
    {
      std::cout << "le elem n'existe pas" << std::endl;
    }

    for(TiXmlElement* sample = elem->FirstChildElement("Class"); sample != NULL; sample = sample->NextSiblingElement())
    {
      // Get the current value of the statistic vector
      int name, value;
      sample->QueryIntAttribute("name", &name);
      sample->QueryIntAttribute("value", &value);
      elmtsInClass[name] = value;
    }
    for(std::map<int, int>::iterator iClass = elmtsInClass.begin(); iClass != elmtsInClass.end(); ++iClass)
    {
      std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    }
    elem = elem->NextSiblingElement();
    
    for(TiXmlElement* sample = elem->FirstChildElement("Id"); sample != NULL; sample = sample->NextSiblingElement())
    {
      // Get the current value of the statistic vector
      double name, value;
      sample->QueryDoubleAttribute("name", &name);
      sample->QueryDoubleAttribute("value", &value);
      polygon[name] = value;
      randomPositionInPolygon[name] = static_cast<int>(generator->GetUniformVariate(0, polygon[name]));
    }
       
    //Initialisation counter of pixel raised in each classes
    std::map<int, int> nbPixelsRaised;
    std::map<int, int> nbSamples;
    for(std::map<int, int>::iterator iClass = elmtsInClass.begin(); iClass != elmtsInClass.end(); ++iClass)
    {
      nbSamples[(*iClass).first] = GetParameterInt("samples");
      if(nbSamples[(*iClass).first] > (*iClass).second)
      {
        nbSamples[(*iClass).first] = (*iClass).second;
      }      
    }
    
    /*for(std::map<unsigned long, int>::iterator ipolygon = polygon.begin(); ipolygon != polygon.end(); ++ipolygon)
    {
      std::cout << "Dans le polygon " << (*ipolygon).first << " il y a " << (*ipolygon).second << " pixels." << std::endl;
    }*/
      
    // *** *** 2nd run : SAMPLING   *** ***
    //Loop across tiles
    for(unsigned int row = 0; row < nbTilesY; ++row)
    {
      for(unsigned int column = 0; column < nbTilesX; ++column)
      {         
        //Progression bar printing
        currentProgression = ((row)*(nbTilesY)+(column+1))*10/(nbTilesX*nbTilesY);        
        if(currentProgression > stepsProgression)
        {
          std::cout<<stepsProgression*10<<"%..."<<std::flush;
          stepsProgression++;
        }
          
        //Tiles dimensions
        unsigned long startX = column*sizeTilesX;
        unsigned long startY = row*sizeTilesY;
        unsigned long sizeX  = vcl_min(sizeTilesX,sizeImageX-startX);
        unsigned long sizeY  = vcl_min(sizeTilesY,sizeImageY-startY);
         
        //Extraction of the image
        ExtractROIFilterType::Pointer extractROIFilter = ExtractROIFilterType::New();
        extractROIFilter->SetInput(image);
        extractROIFilter->SetStartX(startX);
        extractROIFilter->SetStartY(startY);
        extractROIFilter->SetSizeX(sizeX);
        extractROIFilter->SetSizeY(sizeY);
        extractROIFilter->Update();           
                
        //Extraction of the shape file
        OGRLinearRing spatialFilterRing;
        OGRPoint ul,ur,lr,ll;
        ImageType::IndexType urIndex;
        ImageType::IndexType llIndex;
        
        urIndex[0] = startX + sizeX;
        urIndex[1] = startY + sizeY;
        llIndex[0] = startX;
        llIndex[1] = startY;        
          
        itk::Point<double, 2> ulPoint;
        itk::Point<double, 2> lrPoint;
          
        image->TransformIndexToPhysicalPoint(urIndex, lrPoint);
        image->TransformIndexToPhysicalPoint(llIndex, ulPoint);  
           
        ul.setX(ulPoint[0]);
        ul.setY(ulPoint[1]);
        lr.setX(lrPoint[0]);
        lr.setY(lrPoint[1]);
          
        otb::ogr::DataSource::Pointer vectorData2 = otb::ogr::DataSource::New(GetParameterString("shp").c_str(), otb::ogr::DataSource::Modes::Read);
        otb::ogr::Layer filtered2 = vectorData2->GetLayer(0);
        filtered2.SetSpatialFilterRect(ul.getX(), ul.getY(), lr.getX(), lr.getY());
        otb::ogr::Layer::const_iterator featIt = filtered2.begin(); 
          
        unsigned int c = 0;
          
        //Loop across the features in the layer
        for(; featIt!=filtered2.end(); ++featIt,++c)
        {                     
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          bool testPoly = false;
          bool testLineBuffers = false;
          if(geom->getGeometryType() == wkbLineString)
          {
            // Il faut redéfinir geom
            //geom = geom->Buffer(3);
            //std::cout<< "Coord : " << geom->getCoordinateDimension() << std::endl;
            testLineBuffers = true;
          }
           
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {
            testPoly = true;
            //std::cout<< "Coord : " << featIt->ogr().GetFID() << std::endl;
          }
              
          //Class name recuperation
          int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());
            
          //We are dealing with simple polygons
          if(testPoly || testLineBuffers)
          {           
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);          
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing(); 
              
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
              
            //Compute the number of pixel we need to sample in each polygons
            int nbPixelsInPolygon = static_cast<int>((nbSamples[className])*(polygon[featIt->ogr().GetFID()])/(elmtsInClass[className]));    
              
            //If this number is less then 1, we force it to be 1
            if(nbPixelsInPolygon < 1)
            {
              nbPixelsInPolygon = 1;
              polyForced++;
            }
              
            //Compute of the period of sampling, every n-pixels we raise one
            int periodOfSampling = static_cast<int>(polygon[featIt->ogr().GetFID()]/nbPixelsInPolygon);    
                          
            //Counters initialisations            
            //Shifted counter, to anticipate the position of the next raised pixel            
            counterPixelsInPolygonShifted[featIt->ogr().GetFID()] = counterPixelsInPolygon[featIt->ogr().GetFID()] + periodOfSampling;
            //Position of the next pixel raised
            int nextPixelRaisedPosition = periodOfSampling;
                                                
            //Loop across pixels in the tile
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {   
              itk::Point<double, 2> point;
              
              //Boolean variable, we sample the current pixel or not
              bool resultTest = false;
                                        
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(), point); 
                
              //Access to the pixel value
              ImageType::PixelType pixelValue = it.Get();
                
              //Tranformation in OGRPoint in order to know if it's in the geometry
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]); 
                
              //Test if the pixel is not a "No-Data-Pixel", with one of its elmts = 0
              bool noDataTest = false;            
              for (unsigned int i=0; i<nbComponents; i++)
              {   
                if(pixelValue[i] == noDataValue)
                {
                  noDataTest = true; 
                }  
              }
                
              //Pre-test is the polygon at least one pixel
              if(polygon[featIt->ogr().GetFID()]!=0)
              {
                //Succession of tests to know in whoch mode we are
                //Exhautive mode : we extract all pixels in every polygons in every classes
                if (samplingMode == "exhaustive")
                {
                  resultTest= true;
                }
                  
                //Random mode : we extract nbSamples pixels randomly in all the image
                if (samplingMode == "random")
                {
                  //The probability of sampling a pixel is function of the number of pixels in every classes
                  float probability = static_cast<float>(nbSamples[className])/static_cast<float>(nbPixelsGlobal);
                  if(generator->GetUniformVariate(0, 1) < probability)
                  {
                    resultTest= true;        
                  }
                }
                  
                //Random mode equally : we extract nbSAmples pixels for each classes
                if (samplingMode == "randomequally")
                {
                  //The probability of sampling a pixel is function of the number of pixels in each classes
                  float probability = static_cast<float>(nbSamples[className])/static_cast<float>(elmtsInClass[className]);
                  if(generator->GetUniformVariate(0, 1) < probability)
                  {
                    resultTest= true;        
                  }
                }
                  
                //Periodic : we extract, more or less, nbsamples pixels for each classes every n pixels
                if (samplingMode == "periodic")
                {
                  //In a polygon where we only need one pixel, we raise it at a radom position
                  if((counterPixelsInPolygon[featIt->ogr().GetFID()] == randomPositionInPolygon[featIt->ogr().GetFID()])&&(nbPixelsInPolygon == 1))
                  {
                    resultTest= true;
                  } 
                  //If we need more then one pixel in the polygon, we sample a pixel periodicly, every n-pixels
                  if((counterPixelsInPolygon[featIt->ogr().GetFID()]%periodOfSampling)==0 && (nbPixelsInPolygon != 1))
                  {
                    resultTest= true;                  
                  }
                }
                  
                //Periodic random : we extract, more or less, nbsmaples pixels for each classes every n+-delta pixels ( 0<delta<n/2 )
                if (samplingMode == "periodicrandom")
                {
                  //In a polygon where we only need one pixel, we raise it at a radom position
                  if((counterPixelsInPolygon[featIt->ogr().GetFID()] == randomPositionInPolygon[featIt->ogr().GetFID()])&&(nbPixelsInPolygon == 1))
                  {
                    resultTest= true;
                  } 
                  //If we need more then one pixel in the polygon
                  else if(nbPixelsInPolygon != 1)
                  {
                    //The first pixel raised is randomly choosen
                    if(counterPixelsInPolygon[featIt->ogr().GetFID()] == static_cast<int>(generator->GetUniformVariate(0, (periodOfSampling/2))))
                    {
                      resultTest= true;  
                    }
                      
                    //We raised the pixel if we are at the good position
                    if(counterPixelsInPolygon[featIt->ogr().GetFID()] == nextPixelRaisedPosition)
                    {
                      resultTest= true; 
                    }
                      
                    //Every n-pixels we compute the position of the next pixel sampled
                    if(counterPixelsInPolygonShifted[featIt->ogr().GetFID()]%periodOfSampling == 0)
                    {
                      int sign = generator->GetUniformVariate(0, 1);
                      int rdm = static_cast<int>(generator->GetUniformVariate(0, (periodOfSampling/2)));   
                      
                      if (sign<0.5)
                      {
                        nextPixelRaisedPosition = counterPixelsInPolygonShifted[featIt->ogr().GetFID()] - rdm;
                      }
                      else
                      {
                        nextPixelRaisedPosition = counterPixelsInPolygonShifted[featIt->ogr().GetFID()] + rdm;
                      }
                    }                                  
                  }
                }
              }
                
              //Test if the current pixel is in a polygon hole
              bool isNotInHole = true;
              
              for (int i=0; i != inPolygon->getNumInteriorRings(); ++i)
              {
                if(inPolygon->getInteriorRing(i) != NULL)
                {
                  OGRLinearRing * interiorRing = inPolygon->getInteriorRing(i);
                  if(interiorRing->isPointInRing(&pointOGR, TRUE))
                  {
                    isNotInHole = false;
                  }
                }
              }  
                
              //If the pixel is not "No-Data" and is in the geometry, them we count it              
              /*&& !(exteriorRing->isPointOnRingBoundary(&pointOGR, TRUE))*/
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE) && isNotInHole )
              {              
                //Test if the current pixel is good to sample or not
                if(resultTest)
                {
                  std::string message = to_string(className);
                  otb::ogr::Feature featureOutput(layer.GetLayerDefn());       
                  
                  //Adding the raised pixel to our output shape file
                  featureOutput.SetGeometry(&pointOGR);     
                    
                  //Completing the text and shape output files with the pixel values
                  for (unsigned int i=0; i<nbComponents; i++)
                  {
                    message += " " + to_string(i+1) + ": " + to_string(pixelValue[i]);  
                    featureOutput.ogr().SetField(i, pixelValue[i]);    
                  }
                    
                  //We also add informations about where the pixel is extract from
                  for(int c = 0; c < preFeature.ogr().GetFieldCount(); ++c)
                  {  
                    if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTString )
                    {
                      featureOutput.ogr().SetField(nbComponents + c, featIt->ogr().GetFieldAsString(c));
                    }  
                    else if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTInteger )
                    {
                      featureOutput.ogr().SetField(nbComponents + c, featIt->ogr().GetFieldAsInteger(c));
                    }
                  }  
                    
                  layer.CreateFeature(featureOutput);
                  myfile << message << std::endl; 
                  //Incrementation of the counter of raised pixels in each classes
                  nbPixelsRaised[className]++;
                }    
                //Incrementation of counters of pixels studied
                counterPixelsInPolygon[featIt->ogr().GetFID()]++;
                counterPixelsInPolygonShifted[featIt->ogr().GetFID()]++;
              }                
            }  
          }   
        }
      }
    } 
    
    myfile.close();
    
    //End of progression bar
    std::cout<<"100%"<<std::endl;    
    //std::cout<<"polyForced" << polyForced<<std::endl;
    //Output the number of pixel sampled in each classes.
    for(std::map<int, int>::iterator iClass = nbPixelsRaised.begin(); iClass != nbPixelsRaised.end(); ++iClass)
    {
      otbAppLogINFO(<< "Number of pixels raised in class " << (*iClass).first << " : "<< (*iClass).second<< std::endl);
    }
    
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::SamplingImageList)