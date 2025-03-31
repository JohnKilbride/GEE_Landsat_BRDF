#!/usr/bin/env python3
"""
Script Name: usage_example.py
Author: John Kilbride
Date: 2025-03-30
Description:
    
    This script demonstrates how to apply BRDF (Bidirectional Reflectance Distribution Function) 
    correction to a Landsat 8 Collection 2 surface reflectance image using Google Earth 
    Engine (GEE) with Python. It uses the apply_cfactor_brdf_correction function 
    from the brdf_correction module to normalize surface reflectance values based 
    on view and illumination geometry.
    
"""

import ee
from brdf_correction.cfactor import apply_cfactor_brdf_correction

ee.Initialize(project='newenglandagb')


def main (username: str):
    
    # Load in a Landast 8 Collection 2 image for demonstration
    demo_image = ee.Image("LANDSAT/LC08/C02/T1_L2/LC08_223064_20150818")

    # Rename the OLI bands to align with the Landsat TM/ETM+ scheme
    oli_old_names = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"]
    oli_new_names = ["B1", "B2", "B3", "B4", "B5", "B7"]
    demo_image = demo_image.select(oli_old_names, oli_new_names)
    
    # Scale the SR values to a range of 0-10000
    # This isn't necessary for the logic to work, you can save the values
    # as floats in the original SR range (0 to 1)
    spectral_bands = demo_image.multiply(0.0000275) \
        .add(-0.2) \
        .multiply(10000) \
        .toInt16()
    demo_image = demo_image.addBands(spectral_bands, None, True)
    
    # Run the BRDF correction
    demo_image_brdf = apply_cfactor_brdf_correction(demo_image)
    
    # Export the image
    task = ee.batch.Export.image.toAsset(
        image = demo_image_brdf.toInt16(),
        description = "BRDF-GEE-Python",
        assetId = f'users/{username}/LC08_L1TP_223064_20150818_20200908_02_T1_PYTHON',
        scale = 30, 
        crs = "EPSG:32622",
        crsTransform = [30,0,628785,0,-30,-523785],
        region = demo_image.geometry(),
        maxPixels = 1e13
        )
    task.start()
    

if __name__ == "__main__":
    
    print("Script start...")
    main(
      username = "JohnBKilbride"
     )
    print("\n...script complete.")
    

      