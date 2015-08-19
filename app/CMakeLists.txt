OTB_CREATE_APPLICATION(NAME otbSampling
                       SOURCES otbSampling.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}
                      )

OTB_CREATE_APPLICATION(NAME BandValues
                       SOURCES bandValues.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}    
                      )

OTB_CREATE_APPLICATION(NAME AnalysisImageList
                       SOURCES analysisImageList.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}                       
                      )
                      
OTB_CREATE_APPLICATION(NAME SamplingImageList
                       SOURCES samplingImageList.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}                       
                      ) 
                      
OTB_CREATE_APPLICATION(NAME StrategyImageList
                       SOURCES strategyImageList.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}                       
                      ) 
                     
OTB_CREATE_APPLICATION(NAME TrainRF
                       SOURCES trainRF.cxx
                       LINK_LIBRARIES ${${otb-module}_LIBRARIES}               
                      )