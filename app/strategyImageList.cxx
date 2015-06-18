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
  
class StrategyImageList : public Application
{
public:
  typedef StrategyImageList Self;
  typedef itk::SmartPointer<Self> Pointer; 
  itkNewMacro(Self);
  itkTypeMacro(StrategyImageList, otb::Application);
  
private:
  void DoInit()
  {
    SetName("StrategyImageList");
    SetDescription("This application generate the right sampling strategy for images from .xml files.");    
    
    AddParameter(ParameterType_InputFilenameList, "xml", "XML Analysis File");
    AddParameter(ParameterType_OutputFilename, "out", "Output XML file");
    AddParameter(ParameterType_Int, "samples", "Number of samples per classes");
    
    AddParameter(ParameterType_Choice, "strategy", "Strategy of sampling");
    AddChoice("strategy.equally", "Equal repartion across the images");
    AddChoice("strategy.proportional", "Proportional repartion across the images");
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  { 
    std::map<int, std::map<int, int> > nbSamples;
    std::map<int, std::map<int, int> > elmtsInClass;
    std::map<int, int> elmtsInClassGlobal;
    int imageCount = 0;
    
    //Strategy recuperation from the parameter
    const std::string samplingStrategy = GetParameterString("strategy");
    
    std::vector<std::string> xmlFilenameList = GetParameterStringList("xml");
    
    for (std::vector<std::string>::iterator f = xmlFilenameList.begin() ; f != xmlFilenameList.end(); ++f)
    {
      std::cout << *f << std::endl;
      
      TiXmlDocument doc(*f);
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
        elmtsInClassGlobal[name] += value;
        elmtsInClass[imageCount][name] = value;
      }   
      
      imageCount++;
    }
        
    for(std::map<int, std::map<int, int> >::iterator iImage = elmtsInClass.begin(); iImage != elmtsInClass.end(); ++iImage)
    {
      for(std::map<int, int>::iterator iClass = (*iImage).second.begin(); iClass != (*iImage).second.end(); ++iClass)
      {
        if (samplingStrategy == "equally")
        {
          nbSamples[(*iImage).first][(*iClass).first] = GetParameterInt("samples");
          if(nbSamples[(*iImage).first][(*iClass).first] > (*iClass).second)
          {
            nbSamples[(*iImage).first][(*iClass).first] = (*iClass).second;
          }
        }
        
        else if (samplingStrategy == "proportional")
        {
          
          std::cout << "coef " << static_cast<float>((*iClass).second) / static_cast<float>(elmtsInClassGlobal[(*iClass).first]) << std::endl;
          nbSamples[(*iImage).first][(*iClass).first] = GetParameterInt("samples") * static_cast<float>((*iClass).second) / static_cast<float>(elmtsInClassGlobal[(*iClass).first]);
          
          if(nbSamples[(*iImage).first][(*iClass).first] > (*iClass).second)
          {
            nbSamples[(*iImage).first][(*iClass).first] = (*iClass).second;
          }
        } 
      }        
    }   
    
    TiXmlDocument docGlobal;
    TiXmlDeclaration* decl = new TiXmlDeclaration("1.0", "", "");
    docGlobal.LinkEndChild(decl);
      
    TiXmlElement * root = new TiXmlElement("StrategyGlobal");
    docGlobal.LinkEndChild(root);     
    
    for(std::map<int, std::map<int, int> >::iterator iImage = nbSamples.begin(); iImage != nbSamples.end(); ++iImage)
    {
      TiXmlElement * featureImage = new TiXmlElement("Image");
      featureImage->SetDoubleAttribute("name", (*iImage).first);
      root->LinkEndChild(featureImage); 
      for(std::map<int, int>::iterator iClass = (*iImage).second.begin(); iClass != (*iImage).second.end(); ++iClass)
      {
        TiXmlElement * featureClass = new TiXmlElement("Class");
        featureClass->SetDoubleAttribute("name", (*iClass).first);
        featureClass->SetAttribute("value", (*iClass).second);
        featureImage->LinkEndChild(featureClass);
      }            
    }
    
    docGlobal.SaveFile(GetParameterString("out").c_str());
         
  }    
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::StrategyImageList)