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
#include "otbStreamingTraits.h"
#include "otbStreamingMinMaxVectorImageFilter.h"
#include "otbStandardFilterWatcher.h"
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
  
class BandValues : public Application
{
public:
  typedef BandValues Self;
  typedef itk::SmartPointer<Self> Pointer; 
  
  typedef FloatVectorImageType                 ImageType;  
  typedef ImageType::InternalPixelType         ImagePixelType;
  typedef UInt32ImageType                      LabelImageType;
  typedef LabelImageType::InternalPixelType    LabelImagePixelType; 
  typedef itk::ImageRegionIterator<ImageType>  IteratorType;
  typedef otb::ImageFileReader<ImageType>      ReaderType;
  typedef otb::StreamingMinMaxVectorImageFilter<ImageType> StreamingMinMaxVectorImageFilterType;
  ImageType::IndexType IndexType;
  
  typename ImageType::RegionType            polygonRegion;
  
  typedef otb::MultiChannelExtractROI<ImagePixelType, ImagePixelType>  ExtractROIFilterType;
  
  itkNewMacro(Self);

  itkTypeMacro(BandValues, otb::Application);
  
private:
  void DoInit()
  {
    SetName("BandValues");
    SetDescription("This application output the min and max values for each band of the input image. ");
    
    AddParameter(ParameterType_InputImage, "in", "Input Image");    
    
  }

  void DoUpdateParameters()
  {
  }

  void DoExecute()
  {  
    //Input image
    ImageType::Pointer image = GetParameterImage("in");
    image->UpdateOutputInformation();    
    
    StreamingMinMaxVectorImageFilterType::Pointer filter = StreamingMinMaxVectorImageFilterType::New();
     
    //ReaderType::Pointer reader = ReaderType::New();
    //reader->SetFileName(image);
    
    filter->GetStreamer()->SetNumberOfLinesStrippedStreaming( 10 );
    filter->SetInput(image);
    otb::StandardFilterWatcher watcher(filter, "Min Max Computation");
    filter->Update();
    
    std::cout << filter->GetMaximum()[0] << " " << filter->GetMaximum()[1]<< " " << filter->GetMaximum()[2]<< " " << filter->GetMaximum()[3] << std::endl;
   
   
  }
};
}
}

OTB_APPLICATION_EXPORT(otb::Wrapper::BandValues)