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
#include <iterator>     // std::advance

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
  
class otbSampling : public Application
{
public:
  typedef otbSampling Self;
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

  itkTypeMacro(otbSampling, otb::Application);
  
private:
  void DoInit()
  {
    SetName("otbSampling");
    SetDescription("This application generates a .txt file from an image and a vectoriel file. ");
    
    AddParameter(ParameterType_InputImage, "in", "Input Image");    
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
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  {  
    ImageType::Pointer image = GetParameterImage("in");
    image->UpdateOutputInformation();    
    
    std::ofstream myfile;
    myfile.open (GetParameterString("out").c_str());
    
    otb::ogr::DataSource::Pointer vectorData = otb::ogr::DataSource::New(GetParameterString("shp").c_str(), otb::ogr::DataSource::Modes::Read);
    
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
    
    int noDataValue = GetParameterInt("nd");
    
    int nbSamples = GetParameterInt("samples");
    
    int sizeTiles = GetParameterInt("tiles");
        
    unsigned int nbComponents = image->GetNumberOfComponentsPerPixel();
    
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
    unsigned long sizeImageX = image->GetLargestPossibleRegion().GetSize()[0];
    unsigned long sizeImageY = image->GetLargestPossibleRegion().GetSize()[1];    
    unsigned int nbTilesX    = sizeImageX/sizeTilesX + (sizeImageX%sizeTilesX > 0 ? 1 : 0);
    unsigned int nbTilesY    = sizeImageY/sizeTilesY + (sizeImageY%sizeTilesY > 0 ? 1 : 0);
 
    
    int nbPixelsGlobal = 0; 
        
    std::map<int, int> elmtsInClass;
    std::map<unsigned long, int> polygon;
    
    std::cout << " -*-*-*-   1ère Passe   -*-*-*- " << std::endl;
    
    // *** *** 1ère passe : PROSPECTION      *** ***    
    for(unsigned int row = 0; row < nbTilesY; ++row)
    {
      for(unsigned int column = 0; column < nbTilesX; ++column)
      {      
        unsigned long startX = column*sizeTilesX;
        unsigned long startY = row*sizeTilesY;
        unsigned long sizeX  = vcl_min(sizeTilesX+1,sizeImageX-startX+1);
        unsigned long sizeY  = vcl_min(sizeTilesY+1,sizeImageY-startY+1);
                
        ExtractROIFilterType::Pointer extractROIFilter = ExtractROIFilterType::New();
        extractROIFilter->SetInput(image);
        extractROIFilter->SetStartX(startX);
        extractROIFilter->SetStartY(startY);
        extractROIFilter->SetSizeX(sizeX);
        extractROIFilter->SetSizeY(sizeY);
        extractROIFilter->Update();
                
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
             
        for(; featIt!=filtered.end(); ++featIt)
        {          
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {   
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing ();    
                  
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
                                  
            int nbOfPixelsInGeom = 0;
            
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {                          
              itk::Point<double, 2> point;
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(), point);           
              
              ImageType::PixelType pixelValue = it.Get();
                            
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]);
                            
              bool noDataTest = false;            
              for (unsigned int i=0; i<nbComponents; i++)
              {   
                if(noDataTest && (pixelValue[i] == noDataValue))
                {
                  noDataTest = true; 
                }  
              }              
                              
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE))
              {
                nbOfPixelsInGeom++;
                nbPixelsGlobal++;
              }            
            }            
            int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
            //std::cout<< " Nb of Pixels in " << className << " is : " << nbOfPixelsInGeom << std::endl;
            
            polygon[featIt->ogr().GetFID()] += nbOfPixelsInGeom;
            elmtsInClass[className] = elmtsInClass[className] + nbOfPixelsInGeom;
          }        
        }
      }     
    }
    
    /* TRACES */
    //
    //std::cout<< "Nb de classes : " << elmtsInClass.size() << std::endl;
    
    //std::cout << "Nb nbPixelsGlobal " << nbPixelsGlobal << std::endl;
    
    //for(std::map<int, int>::iterator iClass = elmtsInClass.begin(); iClass != elmtsInClass.end(); ++iClass)
    //{
    //  std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    //}
    
    std::cout<< "Nb de polygones : " << polygon.size() << std::endl;

    //for(std::map<unsigned long, int>::iterator ipolygon = polygon.begin(); ipolygon != polygon.end(); ++ipolygon)
    //{
    //  std::cout << "Dans le polygon " << (*ipolygon).first << " il y a " << (*ipolygon).second << " pixels." << std::endl;
    //} 
    
    
    std::cout << " -*-*-*-   2ème Passe   -*-*-*- " << std::endl;
    
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize();
    
    const std::string samplingMode = GetParameterString("mode");
             
    int nbPixelsRaised[elmtsInClass.size()+1];
    for(int b=0; b < elmtsInClass.size()+1; b++)
    {
      nbPixelsRaised[b]=0;
    }
    
    // *** *** 2ème passe : ECHANTILLONAGE   *** ***
    for(unsigned int row = 0; row < nbTilesY; ++row)
    {
      for(unsigned int column = 0; column < nbTilesX; ++column)
      {      
        unsigned long startX = column*sizeTilesX;
        unsigned long startY = row*sizeTilesY;
        unsigned long sizeX  = vcl_min(sizeTilesX,sizeImageX-startX);
        unsigned long sizeY  = vcl_min(sizeTilesY,sizeImageY-startY);
                
        ExtractROIFilterType::Pointer extractROIFilter = ExtractROIFilterType::New();
        extractROIFilter->SetInput(image);
        extractROIFilter->SetStartX(startX);
        extractROIFilter->SetStartY(startY);
        extractROIFilter->SetSizeX(sizeX);
        extractROIFilter->SetSizeY(sizeY);
        extractROIFilter->Update();                

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
                
        for(; featIt!=filtered.end(); ++featIt)
        {                     
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          
          int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
               
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);          
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing ();       
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
                                             
            int nbPixelsInPolygon = int ((nbSamples)*(polygon[featIt->ogr().GetFID()])/(elmtsInClass[className]));            
            if(nbPixelsInPolygon < 1)
            {
              nbPixelsInPolygon = 1;
            }
            
            int periodOfSampling = int(polygon[featIt->ogr().GetFID()]/nbPixelsInPolygon);       
            int randomPositionInPolygon = int(generator->GetUniformVariate(0, polygon[featIt->ogr().GetFID()]));
            
            int counterPixelsInPolygon = 0;
            int counterPixelsInPolygonShited = periodOfSampling;
            int nextPixelRaisedPosition = periodOfSampling;
            
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {   
              itk::Point<double, 2> point;
                            
              bool resultTest = false;
                                         
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(), point); 
              ImageType::PixelType pixelValue = it.Get();
                      
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]); 
                             
              bool noDataTest = false;            
              for (unsigned int i=0; i<nbComponents; i++)
              {   
                if(pixelValue[i] == noDataValue)
                {
                  noDataTest = true; 
                }  
              }
                  
              if (samplingMode == "exhaustive")
              {
                resultTest= true;
              }
              
              if (samplingMode == "random")
              {
                float proba = float(nbSamples)/float(nbPixelsGlobal);
                if(generator->GetUniformVariate(0, 1) < proba)
                {
                  resultTest= true;        
                }
              }
              
              if (samplingMode == "randomequally")
              {
                float proba = (float(nbSamples))/(float(elmtsInClass[className]));
                if(generator->GetUniformVariate(0, 1) < proba)
                {
                  resultTest= true;        
                }
              }
              
              if (samplingMode == "periodic")
              {
                if(polygon[featIt->ogr().GetFID()]!=0)
                {
                  if((counterPixelsInPolygon == polygon[featIt->ogr().GetFID()]/2)&&(nbPixelsInPolygon == 1))
                  {
                    resultTest= true;
                  }              
                  else if((counterPixelsInPolygon%periodOfSampling)==0 && (nbPixelsInPolygon != 1))
                  {
                    resultTest= true;                  
                  }
                }
              }
              
              if (samplingMode == "periodicrandom")
              {
                if(polygon[featIt->ogr().GetFID()]!=0)
                {
                  if((counterPixelsInPolygon == randomPositionInPolygon) && (nbPixelsInPolygon == 1))
                  {
                    resultTest= true;
                  }              
                  else if(nbPixelsInPolygon != 1)
                  {
                    if(counterPixelsInPolygon== int(generator->GetUniformVariate(0, (periodOfSampling/2))))
                    {
                      resultTest= true;  
                    }
                    
                    if(counterPixelsInPolygon == nextPixelRaisedPosition)
                    {
                      resultTest= true; 
                    }
                    
                    if(counterPixelsInPolygonShited%periodOfSampling == 0)
                    {
                      int sign = generator->GetUniformVariate(0, 1);
                      int rdm = int(generator->GetUniformVariate(0, (periodOfSampling/2)));   
                      
                      if (sign<0.5)
                      {
                        nextPixelRaisedPosition = counterPixelsInPolygonShited - rdm;
                      }
                      else
                      {
                        nextPixelRaisedPosition = counterPixelsInPolygonShited + rdm;
                      }
                    }                                  
                  }
                }
              }
              
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE))
              {   
                if(resultTest)
                {
                  std::string message = to_string(className);
                  otb::ogr::Feature dstFeature(layer.GetLayerDefn());       
                  
                  dstFeature.SetGeometry(&pointOGR);     
                            
                  for (unsigned int i=0; i<nbComponents; i++)
                  {
                    message += " " + to_string(i+1) + ": " + to_string(pixelValue[i]);  
                    dstFeature.ogr().SetField(i, pixelValue[i]);    
                  }
                  
                  for(int c = 0; c < preFeature.ogr().GetFieldCount(); ++c)
                  {  
                    if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTString )
                    {
                      dstFeature.ogr().SetField(nbComponents + c, featIt->ogr().GetFieldAsString(c));
                    }  
                    else if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTInteger )
                    {
                      dstFeature.ogr().SetField(nbComponents + c, featIt->ogr().GetFieldAsInteger(c));
                    }
                  }  
                  
                  layer.CreateFeature(dstFeature);
                  myfile << message << std::endl;  
                  nbPixelsRaised[className]++;
                }                
                counterPixelsInPolygon++; 
                counterPixelsInPolygonShited++;
              }                
            }  
          }   
        }
      }
    }      
    myfile.close();
    
    for(int b=1; b < elmtsInClass.size()+1; b++)
    {
      std::cout<< "Number of pixels raised : "<< nbPixelsRaised[b]<< std::endl;
    }    
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::otbSampling)