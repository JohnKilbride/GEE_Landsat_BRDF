#!/usr/bin/env python3
"""
Script Name: cfactor_brdf.py
Author: John Kilbride
Date: 2025-03-30
Description:
    
    This function implements a Python version of the c-factor BRDF logic for
    usage in Google Earth Engine. This logic is designed to be mapped over a 
    collection of Landsat images.

    This Python code was developed based on a Google Earth Engine JavaScript
    implementation of the c-factor BRDF correction. 
    
    Original JavaScript: https://code.earthengine.google.com/3a6761dea6f1bf54b03de1b84dc375c6
    Original Authors: 
    
"""

import ee
import math

from typing import Union

def apply_cfactor_brdf_correction(image: ee.Image) -> ee.Image:
    """
    Applies a c-factor BRDF correction to an input Landsat reflectance image, generating
    nadir BRDF-adjusted reflectance (NBAR). This function follows the approach described
    in Roy et al. (2016), using a fixed set of MODIS-derived BRDF spectral model parameters
    (RossThick and LiThin) to normalize each band to nadir view and a specified solar zenith.

    The function:
      1. Computes image geometry (view angles) and solar positioning.
      2. Calculates RossThick and LiThin BRDF parameters for observed and nadir/solar geometry.
      3. Derives a correction factor for each band based on these parameters.
      4. Returns the corrected image in its original band order.
      
     Usage notes:
         1. This function assumes the Landsat data have a Thematic Mapper/Enhanced
            Thematic Mapper + naming scheme. If you are using this logic with the
            OLI/OLI2 sensors you need to rename the spectral bands to align with the 
            B1, B2, B3, B4, B5, and B7 naming convention.
    Args:
        image (ee.Image): The input Landsat reflectance image to be corrected.

    Returns:
        ee.Image: A BRDF-corrected Landsat reflectance image (NBAR) with the same band order.
    """
    # Constants & coefficients
    coeffs_by_band = {
        'B1': {'fiso': 0.0774, 'fgeo': 0.0079, 'fvol': 0.0372},
        'B2': {'fiso': 0.1306, 'fgeo': 0.0178, 'fvol': 0.0580},
        'B3': {'fiso': 0.1690, 'fgeo': 0.0227, 'fvol': 0.0574},
        'B4': {'fiso': 0.3093, 'fgeo': 0.0330, 'fvol': 0.1535},
        'B5': {'fiso': 0.3430, 'fgeo': 0.0453, 'fvol': 0.1154},
        'B7': {'fiso': 0.2658, 'fgeo': 0.0387, 'fvol': 0.0639}
    }

    # Collect original band order so we can return them correctly at the end
    input_band_names = image.bandNames()

    # Compute geometry, view angles, and solar 
    corners = find_corners(image)
    image = view_angles(image, corners)
    image = solar_position(image)
    image = sun_zen_out(image)


    # Compute RossThick and LiThin parameters
    image = image.addBands(
        image.expression(
            "i.sunAz - i.viewAz", 
            {'i': image}
        ).rename('relativeSunViewAz'),
        overwrite=True
    )
        
    # RossThick
    image = ross_thick(image, 'kvol', 'sunZen', 'viewZen', 'relativeSunViewAz')
    image = ross_thick(image, 'kvol0', 'sunZenOut', 0, 0)

    # LiThin
    image = li_thin(image, 'kgeo', 'sunZen', 'viewZen', 'relativeSunViewAz')    
    image = li_thin(image, 'kgeo0', 'sunZenOut', 0, 0)
    
    # Apply the adjustment to each band
    for band_name, cdict in coeffs_by_band.items():
        image = apply_c_factor(image, band_name, cdict)

    return image.select(input_band_names)


def view_angles(
    image: ee.Image, 
    corners: dict
    ) -> ee.Image:
    """
    Compute satellite viewing azimuth and zenith angles based on image footprint.

    Adds:
        - viewAz: viewing azimuth angle in radians, orthogonal to footprint slope
        - viewZen: viewing zenith angle in radians, estimated from pixel distance to image edges
    """
    # Hard coded parameters
    max_dist = 1000000
    max_sat_zen = 7.5

    # Compute the center points of the upper and lower footprint corners
    upper_center = point_between(corners['upperLeft'], corners['upperRight'])
    lower_center = point_between(corners['lowerLeft'], corners['lowerRight'])

    # Compute footprint slope
    slope_val = compute_slope(lower_center, upper_center)
    slope_perp = ee.Number(-1).divide(slope_val)    

    # Add viewAz
    image = image.addBands(
        image.expression(
            "PI/2 - atan(slope_perp)",
            {'slope_perp': slope_perp, 'PI': math.pi}
        ).rename('viewAz'),
        overwrite=True
    )

    # Compute the center points of the Left and Right hand sides of the foot print
    left_line  = to_line(corners['upperLeft'], corners['lowerLeft'])
    right_line = to_line(corners['upperRight'], corners['lowerRight'])

    # Create distance images from the lines
    left_dist  = ee.FeatureCollection(left_line).distance(max_dist)
    right_dist = ee.FeatureCollection(right_line).distance(max_dist)

    # Compute the view angle zenith image
    view_zen = right_dist.multiply(max_sat_zen*2) \
        .divide(right_dist.add(left_dist)) \
        .subtract(max_sat_zen)

    # Convert to radians
    image = image.addBands(
        view_zen.multiply(math.pi).divide(180).rename('viewZen'),
        overwrite = True
    )
    
    return image


def solar_position(
    image: ee.Image
    ) -> ee.Image:
    """
    Compute solar zenith and azimuth angles based on image time, location, and date.

    Adds the following bands:
        - latRad, longDeg: latitude in radians, longitude in degrees
        - hourGMT, jdp, jdpr: time in hours GMT, Julian date proportion, and radians
        - meanSolarTime, localSolarDiff, trueSolarTime: time-related solar quantities
        - angleHour, delta: hour angle and solar declination
        - cosSunZen, sunZen: cosine and angle of solar zenith
        - sinSunAzSW, cosSunAzSW, sunAzSW: azimuth components from south-west reference
        - sunAz: final solar azimuth angle in radians [0, 2π]
    """
    # Convert the start time to date
    date = ee.Date(ee.Number(image.get('system:time_start')))
    sec_in_hr = 3600

    # Add longitude and latitude in degrees/radians
    lonlat = ee.Image.pixelLonLat()
    image = image.addBands(
        lonlat.select('longitude').rename('longDeg'), overwrite=True
    )
    image = image.addBands(
        lonlat.select('latitude')
        .multiply(math.pi).divide(180)
        .rename('latRad'),
        overwrite=True
    )

    # Hour in GMT and Julian Date Proportion
    hour_gmt = date.getRelative('second', 'day').divide(sec_in_hr)
    jdp = date.getFraction('year')
    image = image.addBands(ee.Image.constant(hour_gmt).rename('hourGMT'), overwrite=True)
    image = image.addBands(ee.Image.constant(jdp).rename('jdp'), overwrite=True)

    # Julian Date Proportion in radians
    image = image.addBands(
        image.expression(
            'jdp * 2 * PI',
            {'jdp': image.select('jdp'), 'PI': math.pi}
        ).rename('jdpr'),
        overwrite=True
    )

    # meanSolarTime
    image = image.addBands(
        image.expression(
            'hourGMT + longDeg / 15',
            {
                'hourGMT': image.select('hourGMT'),
                'longDeg': image.select('longDeg')
            }
        ).rename('meanSolarTime'),
        overwrite=True
    )

    # localSolarDiff
    image = image.addBands(
        image.expression(
            '(0.000075 + 0.001868 * cos(jdpr) - 0.032077 * sin(jdpr)'
            '+ -0.014615 * cos(2 * jdpr) - 0.040849 * sin(2 * jdpr))'
            '* 12 * 60 / PI',
            {
                'jdpr': image.select('jdpr'),
                'PI': math.pi
            }
        ).rename('localSolarDiff'),
        overwrite=True
    )

    # trueSolarTime
    image = image.addBands(
        image.expression(
            'meanSolarTime + localSolarDiff / 60 - 12',
            {
                'meanSolarTime': image.select('meanSolarTime'),
                'localSolarDiff': image.select('localSolarDiff')
            }
        ).rename('trueSolarTime'),
        overwrite=True
    )

    # angleHour
    image = image.addBands(
        image.expression(
            'trueSolarTime * 15 * PI / 180',
            {
                'trueSolarTime': image.select('trueSolarTime'),
                'PI': math.pi
            }
        ).rename('angleHour'),
        overwrite=True
    )

    # Solar declination (delta)
    image = image.addBands(
        image.expression(
            '0.006918 - 0.399912 * cos(jdpr) + 0.070257 * sin(jdpr)'
            ' - 0.006758 * cos(2 * jdpr) + 0.000907 * sin(2 * jdpr)'
            ' - 0.002697 * cos(3 * jdpr) + 0.001480 * sin(3 * jdpr)',
            {'jdpr': image.select('jdpr')}
        ).rename('delta'),
        overwrite=True
    )

    # cosSunZen
    image = image.addBands(
        image.expression(
            'sin(latRad) * sin(delta) + cos(latRad) * cos(delta) * cos(angleHour)',
            {
                'latRad': image.select('latRad'),
                'delta': image.select('delta'),
                'angleHour': image.select('angleHour')
            }
        ).rename('cosSunZen'),
        overwrite=True
    )

    # sunZen
    image = image.addBands(
        image.expression(
            'acos(cosSunZen)',
            {
                'cosSunZen': image.select('cosSunZen')
            }
        ).rename('sunZen'),
        overwrite=True
    )

    # sinSunAzSW (clamped)
    sin_sun_az_sw = image.expression(
        'cos(delta) * sin(angleHour) / sin(sunZen)',
        {
            'delta': image.select('delta'),
            'angleHour': image.select('angleHour'),
            'sunZen': image.select('sunZen')
        }
    ).rename('sinSunAzSW')
    sin_sun_az_sw_clamped = sin_sun_az_sw.clamp(-1, 1).rename('sinSunAzSW')

    image = image.addBands(sin_sun_az_sw_clamped, overwrite=True)

    # cosSunAzSW
    image = image.addBands(
        image.expression(
            '(-cos(latRad) * sin(delta) + sin(latRad) * cos(delta) * cos(angleHour))'
            '/ sin(sunZen)',
            {
                'latRad': image.select('latRad'),
                'delta': image.select('delta'),
                'angleHour': image.select('angleHour'),
                'sunZen': image.select('sunZen')
            }
        ).rename('cosSunAzSW'),
        overwrite=True
    )

    # sunAzSW
    image = image.addBands(
        image.expression(
            'asin(sinSunAzSW)',
            {'sinSunAzSW': image.select('sinSunAzSW')}
        ).rename('sunAzSW'),
        overwrite=True
    )

    # Adjust sunAzSW based on cosSunAzSW and sinSunAzSW
    # 1) If cosSunAzSW <= 0 => sunAzSW = PI - sunAzSW
    image = image.addBands(
        image.select('sunAzSW').where(
            image.select('cosSunAzSW').lte(0),
            ee.Image.constant(math.pi).subtract(image.select('sunAzSW'))
        ).rename('sunAzSW'),
        overwrite=True
    )

    # 2) If cosSunAzSW > 0 and sinSunAzSW <= 0 => sunAzSW = 2*PI + sunAzSW
    image = image.addBands(
        image.select('sunAzSW').where(
            image.select('cosSunAzSW').gt(0).And(image.select('sinSunAzSW').lte(0)),
            ee.Image.constant(2 * math.pi).add(image.select('sunAzSW'))
        ).rename('sunAzSW'),
        overwrite=True
    )

    # sunAz = sunAzSW + PI
    image = image.addBands(
        image.expression(
            'sunAzSW + PI',
            {
                'sunAzSW': image.select('sunAzSW'),
                'PI': math.pi
            }
        ).rename('sunAz'),
        overwrite=True
    )

    # If sunAz > 2*PI => sunAz = sunAz - 2*PI
    image = image.addBands(
        image.select('sunAz').where(
            image.select('sunAz').gt(2 * math.pi),
            image.select('sunAz').subtract(2 * math.pi)
        ).rename('sunAz'),
        overwrite=True
    )

    return image


def sun_zen_out(
    image: ee.Image
    ) -> ee.Image:
    """
    Estimate sun zenith angle at the image center ('sunZenOut') 
    using a latitude-based polynomial.

    Uses the centroid latitude of the image footprint and applies an empirical 
    polynomial approximation from NASA references to compute sun zenith 
    angle in radians.
    """
    center_lat = ee.Number(
        ee.Geometry(image.get('system:footprint'))
            .bounds()
            .centroid(30)
            .coordinates()
            .get(0)
        ).multiply(math.pi).divide(180)

    # Expression for sunZenOut is the same as in JS:
    # (31.0076 - 0.1272*lat + 0.01187*lat^2 + ...) * pi/180
    # For demonstration:
    image = image.addBands(
        ee.Image.constant(
            center_lat.expression(
                "(31.0076 - 0.1272*A + 0.01187*pow(A,2) + 2.40E-05*pow(A,3) "
                "- 9.48E-07*pow(A,4) - 1.95E-09*pow(A,5) + 6.15E-11*pow(A,6)) * PI/180",
                {'A': center_lat, 'PI': math.pi}
            )
        ).rename('sunZenOut'),
        overwrite=True
    )
    
    return image


def ross_thick(
    image: ee.Image,
    band_name: str,
    sun_zen: Union[str, float, int],
    view_zen: Union[str, float, int],
    relative_sun_view_az: Union[str, float, int]
    ) -> ee.Image:
    """
    Compute the Ross-Thick BRDF kernel and add the result as `band_name`.

    Implements:
        1) cosPhaseAngle = cos(sunZen) * cos(viewZen)
                           + sin(sunZen) * sin(viewZen) * cos(relativeAzimuth)
        2) phaseAngle = acos(cosPhaseAngle)
        3) band_name = ((π/2 - phaseAngle) * cosPhaseAngle + sin(phaseAngle)) /
                       (cos(sunZen) + cos(viewZen)) - π/4
    """
    # Convert parameters to ee.Images
    sun_zen_img = get_band_or_constant(image, sun_zen)
    view_zen_img = get_band_or_constant(image, view_zen)
    rel_az_img = get_band_or_constant(image, relative_sun_view_az)

    # 1) cosPhaseAngle
    cos_phase_angle = image.expression(
        'cos(sZ) * cos(vZ) + sin(sZ) * sin(vZ) * cos(rAz)',
        {
            'sZ': sun_zen_img,
            'vZ': view_zen_img,
            'rAz': rel_az_img
        }
    ).rename('cosPhaseAngle')
    image = image.addBands(cos_phase_angle, overwrite=True)

    # 2) phaseAngle
    phase_angle = image.expression(
        'acos(cpa)',
        {
            'cpa': cos_phase_angle
        }
    ).rename('phaseAngle')
    image = image.addBands(phase_angle, overwrite=True)

    # 3) RossThick expression
    # '((PI/2 - phaseAngle)*cosPhaseAngle + sin(phaseAngle)) / (cos(sunZen) + cos(viewZen)) - PI/4'
    ross_expr = image.expression(
        '((PI/2 - pa)*cpa + sin(pa)) / (cos(sZ) + cos(vZ)) - PI/4',
        {
            'pa': phase_angle,
            'cpa': cos_phase_angle,
            'sZ': sun_zen_img,
            'vZ': view_zen_img,
            'PI': math.pi
        }
    ).rename(band_name)
    image = image.addBands(ross_expr, overwrite=True)

    return image


def li_thin(
    image: ee.Image,
    band_name: str,
    sun_zen: Union[str, float, int],
    view_zen: Union[str, float, int],
    relative_sun_view_az: Union[str, float, int],
    h_b: float = 2
    ) -> ee.Image:
    """
    Compute the Li-Thin BRDF kernel and add the result as `band_name`.

    Implements:
        1) sunZenPrime = atan(h/b * tan(sunZen))
           viewZenPrime = atan(h/b * tan(viewZen))
        2) cosPhaseAnglePrime = cos(sunZenPrime) * cos(viewZenPrime)
                                 + sin(sunZenPrime) * sin(viewZenPrime) * cos(relativeAzimuth)
        3) distance = sqrt(tan²(sunZenPrime) + tan²(viewZenPrime)
                           - 2 * tan(sunZenPrime) * tan(viewZenPrime) * cos(relativeAzimuth))
        4) temp = 1/cos(sunZenPrime) + 1/cos(viewZenPrime)
        5) cosT = clamp((h/b) * sqrt(distance² + (tan(sunZenPrime)*tan(viewZenPrime)*sin(relativeAzimuth))²) / temp, -1, 1)
           t = acos(cosT)
        6) overlap = (1/π) * (t - sin(t)*cosT) * temp; overlap = 0 if overlap > 0
        7) band_name = overlap - temp + 0.5 * (1 + cosPhaseAnglePrime) * (1/cos(sunZenPrime)) * (1/cos(viewZenPrime))
    """
    # Convert parameters to ee.Image
    sZ = get_band_or_constant(image, sun_zen)
    vZ = get_band_or_constant(image, view_zen)
    rAz = get_band_or_constant(image, relative_sun_view_az)
    hb = ee.Image.constant(h_b)  # h/b ratio

    # Compute sunZenPrime
    sunZenPrime = image.expression(
        'atan(hb * tan(sZ))',
        {'hb': hb, 'sZ': sZ}
    ).rename('sunZenPrime')
    image = image.addBands(sunZenPrime, overwrite=True)

    # Compute viewZenPrime
    viewZenPrime = image.expression(
        'atan(hb * tan(vZ))',
        {'hb': hb, 'vZ': vZ}
    ).rename('viewZenPrime')
    image = image.addBands(viewZenPrime, overwrite=True)

    # Compute cosPhaseAnglePrime
    cosPhaseAnglePrime = image.expression(
        'cos(sZp) * cos(vZp) + sin(sZp) * sin(vZp) * cos(rAz)',
        {
            'sZp': sunZenPrime,
            'vZp': viewZenPrime,
            'rAz': rAz
        }
    ).rename('cosPhaseAnglePrime')
    image = image.addBands(cosPhaseAnglePrime, overwrite=True)

    # Compute distance
    distance = image.expression(
        'sqrt(pow(tan(sZp), 2) + pow(tan(vZp), 2)'
        ' - 2 * tan(sZp) * tan(vZp) * cos(rAz))',
        {
            'sZp': sunZenPrime,
            'vZp': viewZenPrime,
            'rAz': rAz
        }
    ).rename('distance')
    image = image.addBands(distance, overwrite=True)

    # Compute temp
    temp = image.expression(
        '1/cos(sZp) + 1/cos(vZp)',
        {
            'sZp': sunZenPrime,
            'vZp': viewZenPrime
        }
    ).rename('temp')
    image = image.addBands(temp, overwrite=True)

    # Compute cosT (clamped)
    cosT = image.expression(
        'hb * sqrt(pow(distance, 2) + pow(tan(sZp)*tan(vZp)*sin(rAz), 2)) / temp',
        {
            'hb': hb,
            'distance': distance,
            'sZp': sunZenPrime,
            'vZp': viewZenPrime,
            'rAz': rAz,
            'temp': temp
        }
    ).clamp(-1, 1).rename('cosT')
    image = image.addBands(cosT, overwrite=True)

    # Compute acos
    t_img = image.expression(
        'acos(cosT)',
        {'cosT': cosT}
    ).rename('t')
    image = image.addBands(t_img, overwrite=True)

    # Compute the overlap
    overlap = image.expression(
        '(1/PI) * (t - sin(t)*cosT) * temp',
        {
            'PI': math.pi,
            't': t_img,
            'cosT': cosT,
            'temp': temp
        }
    ).rename('overlap')
    image = image.addBands(overlap, overwrite=True)

    # If overlap > 0 => overlap = 0
    overlap_adjusted = overlap.where(overlap.gt(0), 0).rename('overlap')
    image = image.addBands(overlap_adjusted, overwrite=True)

    # 7) Final Li-Thin expression
    li_expr = image.expression(
        'o - tmp + 0.5 * (1 + cpaPrime) * (1/cos(sZp)) * (1/cos(vZp))',
        {
            'o': overlap_adjusted,
            'tmp': temp,
            'cpaPrime': cosPhaseAnglePrime,
            'sZp': sunZenPrime,
            'vZp': viewZenPrime
        }
    ).rename(band_name)
    image = image.addBands(li_expr, overwrite=True)

    return image


def apply_c_factor(
    image: ee.Image,
    band_name: str,
    coefficients: dict
    ) -> ee.Image:
    """
    Apply c-factor BRDF correction to a band using kernel-based coefficients.

    Implements:
        1) brdf('brdf', 'kvol', 'kgeo', coefficients)
        2) brdf('brdf0', 'kvol0', 'kgeo0', coefficients)
        3) cFactor = brdf0 / brdf
        4) band_name = band_name * cFactor

    The final corrected band overwrites `band_name` in the image.
    """
    # Compute BRDF with actual geometry
    image = compute_brdf(image, 'brdf', 'kvol', 'kgeo', coefficients)

    # Compute BRDF with reference geometry
    image = compute_brdf(image, 'brdf0', 'kvol0', 'kgeo0', coefficients)

    # Compute c-factor as ratio of reference to actual BRDF
    c_factor = image.expression(
        'brdf0 / brdf',
        {
            'brdf0': image.select('brdf0'),
            'brdf': image.select('brdf')
        }
    ).rename('cFactor')

    image = image.addBands(c_factor, overwrite=True)

    # Apply c-factor to the original band
    corrected_band = image.expression(
        'orig * cf',
        {
            'orig': image.select(band_name),
            'cf': c_factor
        }
    ).rename(band_name)

    # Overwrite the old spectral bands with the corrected values
    image = image.addBands(corrected_band, overwrite=True)

    return image


def compute_brdf(
    image: ee.Image,
    band_name: str,
    kvol_band: str,
    kgeo_band: str,
    coefficients: dict
    ) -> ee.Image:
    """
    Compute a BRDF-adjusted reflectance band using provided coefficients and kernel bands.
    
    Implements:
      kvol = 3 * i.kvolBand
      kgeo = i.kgeoBand
      brdf = fiso + fvol*kvol + fgeo*kvol

    The result is added to the image with the name given by `band_name`.
    """
    # Access the coefficient values (assuming each is numeric or convertible into ee.Image.constant)
    fiso = ee.Image.constant(coefficients['fiso'])
    fvol = ee.Image.constant(coefficients['fvol'])
    fgeo = ee.Image.constant(coefficients['fgeo'])

    # Prepare the RossThick and LiThin terms from image (kvol, kgeo).
    # NOTE: The commonly used GEE implementation adds "times 3" factor
    #       However, this is a huristic choice that amplifies the adjustment.
    #       In practice, this seems to work well. 
    kvol = image.select(kvol_band).multiply(3)

    # According to the snippet:
    #   '{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}'
    brdf_expr = image.expression(
        'fiso + fvol * kvol + fgeo * kvol',
        {
            'fiso': fiso,
            'fvol': fvol,
            'fgeo': fgeo,
            'kvol': kvol
        }
    ).rename(band_name)

    return image.addBands(brdf_expr, overwrite=True)


def get_band_or_constant(
    image: ee.Image, 
    param: Union[str, float, int]
    ) -> ee.Image:
    """
    Return a band from `image` if `param` is a string, or a 
    constant image if `param` is numeric.
    """
    # Select band if param is a string
    if isinstance(param, str):
        return image.select(param)

    # Otherwise, create constant image from numeric value
    return ee.Image.constant(param)


def find_corners(image: ee.Image) -> dict:
    """
    Extract the four corners (UL, UR, LR, LL) of the image footprint geometry.

    Uses bounding box extents to locate closest matching vertices from the footprint.
    """
    footprint = ee.Geometry(image.get('system:footprint'))
    bounds = ee.List(footprint.bounds().coordinates().get(0))
    coords = footprint.coordinates()
    
    # Get a list of the X and Y coordinates
    xs = coords.map(get_x_coord)
    ys = coords.map(get_y_coord)
    
    # Get the Landsat scene footprint coordinates based on the bounding box corners
    lower_left  = find_corner(get_x_coord(bounds.get(0)), xs, coords)
    lower_right = find_corner(get_y_coord(bounds.get(1)), ys, coords)
    upper_right = find_corner(get_x_coord(bounds.get(2)), xs, coords)
    upper_left  = find_corner(get_y_coord(bounds.get(3)), ys, coords)
    
    return {
        'upperLeft':  upper_left,
        'upperRight': upper_right,
        'lowerRight': lower_right,
        'lowerLeft':  lower_left
        }


def find_corner(
    target_value: ee.Number,
    arr: ee.List,
    coords: ee.List
    ) -> ee.Geometry.Point:
    """
    Return the coordinate from `coords` whose corresponding value in `arr` is
    closest to `target_value`.
    """
    def compute_diff(value):
        return ee.Number(value).subtract(target_value).abs()

    # Compute absolute differences between each value and the target
    diffs = arr.map(compute_diff)

    # Find the minimum difference
    min_val = diffs.reduce(ee.Reducer.min())

    # Get the index of the closest value
    idx = diffs.indexOf(min_val)

    return coords.get(idx)

    
def get_x_coord(pt: ee.Geometry.Point):
    """Extract the x-coordinate from a Point."""
    return ee.Number(ee.List(pt).get(0))


def get_y_coord(pt: ee.Geometry.Point):
    """Extract the y-coordinate from a Point."""
    return ee.Number(ee.List(pt).get(1))


def compute_slope(
    a: ee.Geometry.Point, 
    b: ee.Geometry.Point
    ) -> ee.Number:
    """Compute the slope of the line connecting two Points as dy/dx."""
    y1 = ee.Number(ee.List(a).get(1))
    y2 = ee.Number(ee.List(b).get(1))
    x1 = ee.Number(ee.List(a).get(0))
    x2 = ee.Number(ee.List(b).get(0))
    
    dy = y1.subtract(y2)
    dx = x1.subtract(x2)
    
    return dy.divide(dx)


def point_between(
    pa: ee.Geometry.Point, 
    pb: ee.Geometry.Point
    ) -> ee.Geometry.Point:
    """Return the midpoint between two Points."""
    return ee.Geometry.LineString([pa, pb]).centroid().coordinates()


def to_line(
    pa: ee.Geometry.Point, 
    pb: ee.Geometry.Point
    ) -> ee.Geometry.LineString:
    """Create a LineString from two Points."""
    return ee.Geometry.LineString([pa, pb])


