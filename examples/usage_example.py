import ee

from brdf_correction.cfactor import apply_cfactor_brdf_correction

ee.Initialize(project='newenglandagb')

if __name__ == "__main__":
    
    # Define your GEE username
    username = "JohnBKilbride"
    
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
        crs = demo_image.projection().crs().getInfo(),
        region = demo_image.geometry(),
        maxPixels = 1e13
        )
    task.start()
      