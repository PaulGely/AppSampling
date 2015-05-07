otb_module_test()

otb_test_application(NAME otbSampling
                     APP Convert
                     OPTIONS -in ${INPUTDATA}/QB_Toulouse_Ortho_XS.tif
                             -out ${TEMP}/apTvUtConvertBasicOutput.tif
                     VALID   --compare-image ${NOTOL}
                             ${INPUTDATA}/apTvUtConvertBasicOutput.tif
                       ${TEMP}/apTvUtConvertBasicOutput.tif
)