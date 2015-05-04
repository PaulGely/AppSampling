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
    
    unsigned int nbComp = image->GetNumberOfComponentsPerPixel();
    
    otb::ogr::Layer preFiltered = vectorData->GetLayer(0);    
    otb::ogr::Feature preFeature = preFiltered.GetFeature(0);
    
    for(unsigned int comp = 0; comp<nbComp; ++comp)
    {  
      std::ostringstream fieldoss;
      fieldoss<<"b"<<comp;
      OGRFieldDefn field(fieldoss.str().c_str(), OFTReal);
      layer.CreateField(field, true);
    }  
    
    for(int comp = 0; comp < preFeature.ogr().GetFieldCount(); ++comp)
    {  
      std::ostringstream fieldoss;
      fieldoss<<comp + nbComp;
      OGRFieldDefn field(preFeature.ogr().GetFieldDefnRef(comp)->GetNameRef(), preFeature.ogr().GetFieldDefnRef(comp)->GetType());
      layer.CreateField(field, true);  
    }  
            
    unsigned long sizeTilesX = 500;
    unsigned long sizeTilesY = 500;
    unsigned long sizeImageX = image->GetLargestPossibleRegion().GetSize()[0];
    unsigned long sizeImageY = image->GetLargestPossibleRegion().GetSize()[1];    
    unsigned int nbTilesX    = sizeImageX/sizeTilesX + (sizeImageX%sizeTilesX > 0 ? 1 : 0);
    unsigned int nbTilesY    = sizeImageY/sizeTilesY + (sizeImageY%sizeTilesY > 0 ? 1 : 0);

    //std::cout << "Number of tiles: " << nbTilesX <<" x "<< nbTilesY << std::endl;    
    
    int test = 0;
    int nbPixelCount = 0; 
        
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
                
        //std::cout << "** 1 **  Tiles nb : " << (row)*(nbTilesY)+(column+1) << " over " << nbTilesY*nbTilesX << std::endl;    

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
  
        //std::cout<<"Feature count: "<<filtered.GetFeatureCount(true)<<std::endl;          
  
        otb::ogr::Layer::const_iterator featIt = filtered.begin(); 
             
        for(; featIt!=filtered.end(); ++featIt)
        {          
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {   
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing ();    
                  
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
                                  
            int nbOfPixels = 0;
            
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {                          
              itk::Point<double, 2> point;
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(), point);           
              
              ImageType::PixelType pixelValue = it.Get();
              
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]);
                            
              bool noDataTest = false;            
              for (unsigned int i=0; i<nbComp; i++)
              {   
                if(noDataTest && (pixelValue[i] == noDataValue))
                {
                  noDataTest = true; 
                }  
              }              
                              
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE))
              {
                nbOfPixels++;
                nbPixelCount++;
              }            
            }            
            int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
            //std::cout<< " Nb of Pixels in " << className << " is : " << nbOfPixels << std::endl;
            
            polygon[featIt->ogr().GetFID()] += nbOfPixels;
            elmtsInClass[className] = elmtsInClass[className] + nbOfPixels;
          }        
        }
      }     
        //std::cout << "*** Flag = " << flag << std::endl;
    }
    
    /* TRACES */
    /*
    std::cout<< "Nb de classes : " << elmtsInClass.size() << std::endl;
    
    std::cout << "Nb nbPixelCount " << nbPixelCount << std::endl;
    
    for(std::map<int, int>::iterator iClass = elmtsInClass.begin(); iClass != elmtsInClass.end(); ++iClass)
    {
      std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    }
    
    std::cout<< "Nb de polygones : " << polygon.size() << std::endl;

    for(std::map<unsigned long, int>::iterator ipolygon = polygon.begin(); ipolygon != polygon.end(); ++ipolygon)
    {
      std::cout << "Dans le polygon " << (*ipolygon).first << " il y a " << (*ipolygon).second << " pixels." << std::endl;
    } 
    */
    
    std::cout << " -*-*-*-   2ème Passe   -*-*-*- " << std::endl;
    
    typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
    GeneratorType::Pointer generator = GeneratorType::New();
    generator->Initialize();
             
    int counterLeves[elmtsInClass.size()+1];
    for(int b=0; b < elmtsInClass.size()+1; b++)
    {
      counterLeves[b]=0;
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
                
        //std::cout << "** 2 **  Tiles nb : " << (row)*(nbTilesY)+(column+1) << "/" << nbTilesX*nbTilesY << std::endl;  
        //std::cout << "* X : " << row << "  Y : " << column << std::endl;

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
  
        //std::cout<<"Feature count: "<<filtered.GetFeatureCount(true)<<std::endl;
        otb::ogr::Layer::const_iterator featIt = filtered.begin(); 
        
        for(; featIt!=filtered.end(); ++featIt)
        {                     
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          //std::cout << "*" << geom->getGeometryName () << std::endl;
          
          int className = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
               
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);          
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing ();       
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
             
            int pixC = 0;
            
            int counter=0;
            
            int nbSamples = 2000;
            
            int nbPixelsInPolygon = int ((nbSamples)*(polygon[featIt->ogr().GetFID()])/(elmtsInClass[className]));            
            if(nbPixelsInPolygon < 1)
            {
              nbPixelsInPolygon = 1;
            }
            
            int N = int(polygon[featIt->ogr().GetFID()]/nbPixelsInPolygon);       
                        
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {   
              itk::Point<double, 2> point;
              IteratorType itBis = it;
                            
              bool resultTest = false;
              
              if(polygon[featIt->ogr().GetFID()]!=0)
              {
                if((counter == polygon[featIt->ogr().GetFID()]/2) && (nbPixelsInPolygon == 1))
                {
                  resultTest= true;
                }              
                else if((counter%N)==0 && (nbPixelsInPolygon != 1))
                {
                  resultTest= true; 
                  //int(generator->GetUniformVariate(0, (N/2)))
                  int sign = generator->GetUniformVariate(0, 1);
                  if (sign<0.5)
                  {
                    for(int i=0; i<1; i++)
                    {
                      --itBis;
                    }
                  }
                  else
                  {
                    for(int i=0; i<1; i++)
                    {
                      ++itBis;
                    }
                  }
                  
                }
              }
              
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(itBis.GetIndex(), point); 
              ImageType::PixelType pixelValue = itBis.Get();
                      
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]); 
                             
              bool noDataTest = false;            
              for (unsigned int i=0; i<nbComp; i++)
              {   
                if(pixelValue[i] == noDataValue)
                {
                  noDataTest = true; 
                }  
              }
                            
              /* Random mode */
              //float proba = float(200)/float(nbPixelCount);
              /* Random mode, equally distributed according to the classes. */
              //float proba = (float(2000))/(float(elmtsInClass[className])/**elmtsInClass.size()*/);
              /*
              if(generator->GetUniformVariate(0, 1) < proba)
              {
                resultTest= true;                
                //std::cout<< "Test RDM : " << rdmTest << std::endl;
              }    
              */   
              
              if(!noDataTest && exteriorRing->isPointInRing(&pointOGR, TRUE))
              {   
                if(resultTest)
                {
                  test++;
                  std::string msg = to_string(className);
                  otb::ogr::Feature dstFeature(layer.GetLayerDefn());       
                  
                  dstFeature.SetGeometry(&pointOGR);     
                            
                  for (unsigned int i=0; i<nbComp; i++)
                  {
                    msg += " " + to_string(i+1) + ": " + to_string(pixelValue[i]);  
                    dstFeature.ogr().SetField(i, pixelValue[i]);    
                  }
                  
                  for(int c = 0; c < preFeature.ogr().GetFieldCount(); ++c)
                  {  
                    if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTString )
                    {
                      dstFeature.ogr().SetField(nbComp + c, featIt->ogr().GetFieldAsString(c));
                    }  
                    else if(featIt->ogr().GetFieldDefnRef(c)->GetType() == OFTInteger )
                    {
                      dstFeature.ogr().SetField(nbComp + c, featIt->ogr().GetFieldAsInteger(c));
                    }
                  }  
                  
                  //msg += " " + to_string(point[0]) + " " + to_string(point[1]);
                  layer.CreateFeature(dstFeature);
                  myfile << msg << std::endl;  
                  pixC++;
                  counterLeves[className]++;
                }                
                counter++;    
              }                
            }  
            //std::cout<<"Pix count elmt : "<<elmtsInClass[className]<< std::endl;
            //std::cout<<"Pix  : "<<pixC<< std::endl;
          }   
        }
        //std::cout << "*** Flag = " << flag << std::endl;
      }
    }      
    std::cout << "Nb test " << test << std::endl;
    myfile.close();
    for(int b=1; b < elmtsInClass.size()+1; b++)
    {
      std::cout<< "Counter Pix leves: "<< counterLeves[b]<< std::endl;
    }
    
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::otbSampling)