set(DOCUMENTATION "Sampling application.")

otb_module(AppSampling
  DEPENDS
    OTBApplicationEngine
    OTBGdalAdapters
    OTBITK
    OTBImageBase
    OTBImageIO

  DESCRIPTION
    "${DOCUMENTATION}"
)
