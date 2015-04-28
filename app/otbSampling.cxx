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
#include <vcl_algorithm.h>
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
    
    AddParameter(ParameterType_Int, "nd", "NoData value");
    SetDefaultParameterInt("nd", 0);
    MandatoryOff("nd");
    
    AddParameter(ParameterType_String, "cfield", "Field of class");
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
      fieldoss<<comp;
      OGRFieldDefn field(fieldoss.str().c_str(), OFTReal);
      layer.CreateField(field, true);
    }  
    
    for(int comp = 0; comp < preFeature.ogr().GetFieldCount(); ++comp)
    {  
      std::ostringstream fieldoss;
      fieldoss<<comp + nbComp;
      OGRFieldDefn field(fieldoss.str().c_str(), preFeature.ogr().GetFieldDefnRef(comp)->GetType());
      layer.CreateField(field, true);  
    }  
            
    unsigned long sizeTilesX = 200;
    unsigned long sizeTilesY = 200;
    unsigned long sizeImageX = image->GetLargestPossibleRegion().GetSize()[0];
    unsigned long sizeImageY = image->GetLargestPossibleRegion().GetSize()[1];    
    unsigned int nbTilesX    = sizeImageX/sizeTilesX + (sizeImageX%sizeTilesX > 0 ? 1 : 0);
    unsigned int nbTilesY    = sizeImageY/sizeTilesY + (sizeImageY%sizeTilesY > 0 ? 1 : 0);

    //std::cout << "Number of tiles: " << nbTilesX <<" x "<< nbTilesY << std::endl;    
    
    //unsigned long regionCount = 0;
    int test = 0;
    int nbOfPixelsInImage = sizeImageX * sizeImageY;    
        
    std::map<int, int> elmtOfClass;
    std::map<int, int> poly;
    
    int nbPoly = 0;
    
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
                
        std::cout << "* Tiles nb : " << (row)*(nbTilesY)+(column+1) << std::endl;    

        //int flag = 0;

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
                      
            //flag = 1;
            //myfile << "------ Tiles nb : " << (row)*(nbTilesY)+(column+1) << "------" << std::endl;
            
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
              }     
                            
              //myfile << msg << std::endl;           
            }            
            
            //bool testPressence = false;
            int nomClass = featIt->ogr().GetFieldAsInteger(GetParameterString("cfield").c_str());             
            std::cout<< " Nb of Pixels in " << nomClass << " is : " << nbOfPixels << std::endl;
            
            poly[nbPoly] = nbOfPixels;
            elmtOfClass[nomClass] = elmtOfClass[nomClass] + nbOfPixels;
            
            nbPoly++;
                        
          }        
        }
      }     
        //std::cout << "*** Flag = " << flag << std::endl;
    }
    
    std::cout<< "Nb de classes : " << elmtOfClass.size() << std::endl;
    
    for(std::map<int, int>::iterator iClass = elmtOfClass.begin(); iClass != elmtOfClass.end(); ++iClass)
    {
      std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    }
    
    std::cout<< "Nb de polygones : " << poly.size() << std::endl;

    for(std::map<int, int>::iterator iPoly = poly.begin(); iPoly != poly.end(); ++iPoly)
    {
      std::cout << "Dans le poly " << (*iPoly).first << " il y a " << (*iPoly).second << " pixels." << std::endl;
    } 
         
    std::cout << " -*-*-*-   2ème Passe   -*-*-*- " << std::endl;
    
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
                
        std::cout << "* Tiles nb : " << (row)*(nbTilesY)+(column+1) << "/" << nbTilesX*nbTilesY << std::endl;  
        //std::cout << "* X : " << row << "  Y : " << column << std::endl;

        //int flag = 0;

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
  
        std::cout<<"Feature count: "<<filtered.GetFeatureCount(true)<<std::endl;
        otb::ogr::Layer::const_iterator featIt = filtered.begin(); 
        
        for(; featIt!=filtered.end(); ++featIt)
        {                     
          OGRGeometry * geom = featIt->ogr().GetGeometryRef();
          //std::cout << "*" << geom->getGeometryName () << std::endl;
          
          if(geom->getGeometryType() == wkbPolygon25D || geom->getGeometryType() == wkbPolygon)
          {
            OGRPolygon * inPolygon = dynamic_cast<OGRPolygon *>(geom);          
            OGRLinearRing * exteriorRing = inPolygon->getExteriorRing ();
            //std::cout<<featIt->ogr().GetFieldAsInteger("class")<<std::endl;       
            IteratorType it(extractROIFilter->GetOutput(), extractROIFilter->GetOutput()->GetLargestPossibleRegion());
                      
            //flag = 1;
            //myfile << "------ Tiles nb : " << (row)*(nbTilesY)+(column+1) << "------" << std::endl;
                          
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {   
              itk::Point<double, 2> point;
              extractROIFilter->GetOutput()->TransformIndexToPhysicalPoint(it.GetIndex(), point);           
                
              ImageType::PixelType pixelValue = it.Get();
                      
              OGRPoint pointOGR;
              pointOGR.setX(point[0]);
              pointOGR.setY(point[1]); 
              //std::cout<<pointOGR.getX()<<std::endl;
                             
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
                test++;
                std::string msg = to_string(featIt->ogr().GetFieldAsInteger("class"));
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
              }         
            }
          }   
        }

        //std::cout << "*** Flag = " << flag << std::endl;
      }
    }  
    
    //std::cout << "Nb test " << test << std::endl;
    myfile.close();
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::otbSampling)