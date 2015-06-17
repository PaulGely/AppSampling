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
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  { 
    std::map<int, int> elmtsInClassGlobal;
        
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
      }      
    }
    
    for(std::map<int, int>::iterator iClass = elmtsInClassGlobal.begin(); iClass != elmtsInClassGlobal.end(); ++iClass)
    {
      std::cout << "Dans la classe " << (*iClass).first << " il y a " << (*iClass).second << " pixels." << std::endl;
    }
    
    TiXmlDocument docGlobal;
    TiXmlDeclaration* decl = new TiXmlDeclaration("1.0", "", "");
    docGlobal.LinkEndChild(decl);
      
    TiXmlElement * root = new TiXmlElement("ImageAnalysisGlobal");
    docGlobal.LinkEndChild(root);
      
    for(std::map<int, int>::iterator iClass = elmtsInClassGlobal.begin(); iClass != elmtsInClassGlobal.end(); ++iClass)
    {
            
      TiXmlElement * featureClass = new TiXmlElement("Class");
      featureClass->SetDoubleAttribute("name", (*iClass).first);
      featureClass->SetAttribute("value", (*iClass).second);
      root->LinkEndChild(featureClass);        
    }
    docGlobal.SaveFile(GetParameterString("out").c_str());
  }    
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::StrategyImageList)