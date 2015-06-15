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
  
class AnalysisImageList : public Application
{
public:
  typedef AnalysisImageList Self;
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

  itkTypeMacro(AnalysisImageList, otb::Application);
  
private:
  void DoInit()
  {
    SetName("AnalysisImageList");
    SetDescription("This application analyse the input image, the output is a .xml file.");
    
    AddParameter(ParameterType_InputImage, "in", "Input Image");    
    
    AddParameter(ParameterType_InputFilename, "shp", "Vectoriel File");    
    AddParameter(ParameterType_OutputFilename, "out", "Output Text");       
    AddParameter(ParameterType_String, "cfield", "Field of class");
        
    AddParameter(ParameterType_Int, "tiles", "Size of square tiles");
    SetDefaultParameterInt("tiles", 200);
    MandatoryOff("tiles");
    
    AddParameter(ParameterType_Int, "nd", "NoData value");
    SetDefaultParameterInt("nd", 0);
    MandatoryOff("nd");          
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  {  
    ImageType::Pointer image = GetParameterImage("in");
    image->UpdateOutputInformation(); 
        
    //Input shape file
    otb::ogr::DataSource::Pointer vectorData = otb::ogr::DataSource::New(GetParameterString("shp").c_str(), otb::ogr::DataSource::Modes::Read);
    
    //Output shape file
    std::vector<std::string> options;
    std::string projRef = image->GetProjectionRef();
    OGRSpatialReference oSRS(projRef.c_str());
        
    //No data value for pixels value
    int noDataValue = GetParameterInt("nd");
        
    //Dimension of a side of the square tiles
    int sizeTiles = GetParameterInt("tiles");
    
    //Number of elements in each pixels
    unsigned int nbComponents = image->GetNumberOfComponentsPerPixel();
    
    //Initialisation of the output shape file
    otb::ogr::Layer preFiltered = vectorData->GetLayer(0);    
    otb::ogr::Feature preFeature = preFiltered.GetFeature(0);
        
    unsigned long sizeTilesX = sizeTiles;
    unsigned long sizeTilesY = sizeTiles;
    unsigned long sizeImageX;
    unsigned long sizeImageY;    
    unsigned long nbTilesX;
    unsigned long nbTilesY;
    
    
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
    
    
    otbAppLogINFO(<< "Computing the number of pixels for each polygons and classes" << std::endl);
    
    int polyForced = 0;
    
    int stepsProgression = 0;
    int currentProgression;      
    
    sizeImageX = image->GetLargestPossibleRegion().GetSize()[0];
    sizeImageY = image->GetLargestPossibleRegion().GetSize()[1];
    nbTilesX = sizeImageX/sizeTilesX + (sizeImageX%sizeTilesX > 0 ? 1 : 0);
    nbTilesY = sizeImageY/sizeTilesY + (sizeImageY%sizeTilesY > 0 ? 1 : 0);      
      
    // *** *** 1st run :  PROSPECTION      *** ***    
    int countTest = 0; 
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
        countTest++;
        //std::cout<<"count Test " << countTest << std::endl;
        countTest++;
          
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

        otb::ogr::Layer filtered = vectorData->GetLayer(0);
        filtered.SetSpatialFilterRect(ul.getX(), ul.getY(), lr.getX(), lr.getY());         
            
        otb::ogr::Layer::const_iterator featIt = filtered.begin(); 
          
        //Loop across the features in the layer
        for(; featIt!=filtered.end(); ++featIt)
        {          
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          
          if(!geom)
          {
            std::cout<<featIt->ogr().GetFID()<<std::endl;
            break;  
          } 
              
          bool testPoly = false;
          bool testLineBuffers = false;
          if(geom->getGeometryType() == wkbLineString)
          {
            geom = geom->Buffer(3);
            //std::cout<< "Coord : " << geom->getCoordinateDimension() << std::endl;
            testLineBuffers = true;
          }
            
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {
            testPoly = true;
          }
                     
          //We are dealing with simple polygons
          if(testPoly || testLineBuffers)
          {   
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing();
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
              
            //Number of pixels in a polygon
            int nbOfPixelsInGeom = 0;
              
            //Loop across pixels in the tile
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {                          
              itk::Point<double, 2> point;
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
                if(noDataTest && (pixelValue[i] == noDataValue))
                {
                  noDataTest = true; 
                }  
              }             
                
              //If the pixel is not "No-Data" and is in the geometry, them we count it
              //isPointOnRingBoundary() is not relevent beacause there is not (or very few) pixel excatly on the line boundary...
              /*&& !(exteriorRing->isPointOnRingBoundary(&pointOGR, TRUE))*/
                              
              //Test if the current pixel is in a polygon hole
              bool isNotInHole = true;
              
              for (int i=0; i < inPolygon->getNumInteriorRings(); ++i)
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
                
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE) && isNotInHole)
              {
                nbOfPixelsInGeom++;
                nbPixelsGlobal++;
              }                    
            }
              
            //Class name recuperation
            int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
              
            //Counters update, number of pixel in each classes and in each polygones 
            polygon[featIt->ogr().GetFID()] += nbOfPixelsInGeom;
              
            //Generation of a random number for the sampling in a polygon where we only need one pixel, it's choosen randomly
            elmtsInClass[className] = elmtsInClass[className] + nbOfPixelsInGeom;  
            //std::cout<<"Test"<<std::endl; 
          }   
        }      
      }     
    }
    std::cout<<std::endl;
    
    
    //End of progression bar
    std::cout<<"100%"<<std::endl;
      
    /* TRACES */
    //
    std::cout<< "Nb de classes : " << elmtsInClass.size() << std::endl;
      
    std::cout << "Nb nbPixelsGlobal " << nbPixelsGlobal << std::endl;
      
    /*for(std::map<int, int>::iterator iClass = elmtsInClass.begin(); iClass != elmtsInClass.end(); ++iClass)
    {
      std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    }*/
    for(std::map<int, int>::iterator iImage = elmtsInClass.begin(); iImage != elmtsInClass.end(); ++iImage)
    {      
      std::cout << "Dans l'image la classe " << (*iImage).first << "  a " << (*iImage).second << " pixels." << std::endl;       
    }
      
    std::cout<< "Nb de polygones : " << polygon.size() << std::endl;
     
    /*for(std::map<unsigned long, int>::iterator ipolygon = polygon.begin(); ipolygon != polygon.end(); ++ipolygon)
    {
      std::cout << "Dans le polygon " << (*ipolygon).first << " il y a " << (*ipolygon).second << " pixels." << std::endl;
    }*/   
      
    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild(decl);
      
    TiXmlElement * root = new TiXmlElement("ImageAnalysis");
    doc.LinkEndChild(root);
      
    TiXmlElement * featureImage = new TiXmlElement("Image");
    root->LinkEndChild(featureImage);  
    for(std::map<int, int>::iterator iImage = elmtsInClass.begin(); iImage != elmtsInClass.end(); ++iImage)
    {
            
      TiXmlElement * featureClass = new TiXmlElement("Class");
      featureClass->SetDoubleAttribute("name", (*iImage).first);
      featureClass->SetAttribute("value", (*iImage).second);
      featureImage->LinkEndChild(featureClass);        
    }     
    
    TiXmlElement * featurePolygon = new TiXmlElement("Polygon");
    root->LinkEndChild(featurePolygon);
    for(std::map<unsigned long, int>::iterator ipolygon = polygon.begin(); ipolygon != polygon.end(); ++ipolygon)
    {         
      TiXmlElement * featureId = new TiXmlElement("Id");
      featureId->SetAttribute("name", (*ipolygon).first);
      featureId->SetDoubleAttribute("value", (*ipolygon).second);
      featurePolygon->LinkEndChild(featureId);   
    }
      
    doc.SaveFile(GetParameterString("out").c_str());
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::AnalysisImageList)