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
#include <otbMachineLearningModel.h>
#include "otbRandomForestsMachineLearningModel.h"
#include "itkVariableLengthVector.h"
#include "otbListSampleGenerator.h"
#include <boost/algorithm/string.hpp>

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
  
class TrainRF : public Application
{
public:
  typedef TrainRF Self;
  typedef itk::SmartPointer<Self> Pointer; 
  
  typedef otb::ListSampleGenerator<FloatVectorImageType, VectorDataType> ListSampleGeneratorType;
  
  typedef otb::RandomForestsMachineLearningModel<float, int> RandomForestType;

  typedef otb::MachineLearningModel<float,short>         MachineLearningModelType;  
  typedef MachineLearningModelType::InputSampleType      InputSampleType;
  typedef MachineLearningModelType::InputListSampleType  InputListSampleType;  
  typedef MachineLearningModelType::TargetSampleType     TargetSampleType;
  typedef MachineLearningModelType::TargetListSampleType TargetListSampleType;
  typedef ListSampleGeneratorType::ListLabelType LabelListSampleType;
  typedef ListSampleGeneratorType::LabelType LabelType;
    
  itkNewMacro(Self);

  itkTypeMacro(TrainRF, otb::Application);
  
private:
  void DoInit()
  {
    SetName("TrainRF");
    SetDescription("Train RF. ");
    
    AddParameter(ParameterType_OutputFilename, "in", "Input Image");    
    AddParameter(ParameterType_OutputFilename, "out", "Output model"); 
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  {   
    InputListSampleType::Pointer samples = InputListSampleType::New();
    LabelListSampleType::Pointer labels = LabelListSampleType::New();   
    
    std::ifstream ifs;
    ifs.open(GetParameterString("in").c_str());

    if(!ifs)
    {
      std::cout<<"Could not read file "<<GetParameterString("in")<<std::endl;
    }
        
    unsigned int nbfeatures = 0;

    while (!ifs.eof())
    {
      std::string line;
      std::getline(ifs, line);
      boost::algorithm::trim(line);
        
      if(nbfeatures == 0)
      {
        nbfeatures = std::count(line.begin(),line.end(),' ');
        std::cout<<"Line "<<line<<std::endl;
        std::cout<<"Found "<<nbfeatures<<" features per samples"<<std::endl;
      }      
      
      if(line.size()>1)
      {
        InputSampleType sample(nbfeatures);
        sample.Fill(0);

        std::string::size_type pos = line.find_first_of(" ", 0);  
        
        // Parse label
        LabelType label;
        label[0] = atoi(line.substr(0, pos).c_str());

        bool endOfLine = false;
        unsigned int id = 0;
        
        while(!endOfLine)
        {
          std::string::size_type nextpos = line.find_first_of(" ", pos+1);

          if(pos == std::string::npos)
          {
            endOfLine = true;
            nextpos = line.size()-1;
          }
          
          else
          {          
            std::string feature = line.substr(pos,nextpos-pos);
            std::string::size_type semicolonpos = feature.find_first_of(":");
            id = atoi(feature.substr(0,semicolonpos).c_str());
            sample[id - 1] = atof(feature.substr(semicolonpos+1,feature.size()-semicolonpos).c_str());
             
            pos = nextpos;
            std::cout << "id: " << id <<std::endl;
            std::cout << "sample: " << sample <<std::endl;
          }           

        }      
        
        samples->SetMeasurementVectorSize(itk::NumericTraits<InputSampleType>::GetLength(sample));
        samples->PushBack(sample);
        labels->PushBack(label);
      }     
    }

    std::cout<<"Retrieved "<<samples->Size()<<" samples"<<std::endl;
    ifs.close();
    
    RandomForestType::Pointer classifier = RandomForestType::New();
    classifier->SetInputListSample(samples);
    classifier->SetTargetListSample(labels);
    classifier->SetMaxDepth(5);
    classifier->SetMinSampleCount(10);
    classifier->SetRegressionAccuracy(0.);
    classifier->SetMaxNumberOfCategories(10);
    classifier->SetMaxNumberOfVariables(0);
    classifier->SetMaxNumberOfTrees(100);
    classifier->SetForestAccuracy(0.01);

    classifier->Train();
    classifier->Save(GetParameterString("out")); 
        
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::TrainRF)




