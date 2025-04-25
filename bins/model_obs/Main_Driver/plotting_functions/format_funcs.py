import numpy as np
import math
import re

def format_lat_lon(lat,lon):
   return format_lat_lon_alt(lat,lon,None)

def format_lat_lon_alt(lat,lon,alt):
   latStr = format_lat_or_lon(lat)
   lonStr = format_lat_or_lon(lon,isLatitude=False)
   llaStr = latStr+" "+lonStr
   if (alt):
     altStr = format_number(alt,1) + 'm'
     llaStr = llaStr + " " +altStr
   return llaStr  

def format_number(val,precision):
   return f'{val:.{precision}f}'

"""
 * Format a latitude or longitude value to the given format.  Formats use
 * DD for degrees, MM for minutes, SS for seconds and d, m, s for decimal
 * fractions of degrees, minutes, seconds.  H designates the hemisphere
 * (N,S,E,W).
 * Examples for value -34.496 degrees
 *
 *     DD:MM:SS      ===>  -34:29:45
 *       (if longitude and use360 ===> 326:29:45)
 *     DDH           ===>   34W     (or 34S if longitude)
 *     DD.d          ===>  -34.5
 *     DD.dddH       ===>   34.496W (or 34.496S if longitude)
 *     DD MM" SS.s'  ===>  -34 29" 45.6'
 *
 * value  the value to format
 * format the format
 * isLatitude  true if latitude, false if longitude
 * use360      if true use 0-360 notation instead of -180 to 180 notation
 *
 * @return formatted value
 """
def format_lat_or_lon(value,format='DD.ddH',isLatitude=True,use360=False):
    if (np.isnan(value)):
        return np.nan;
    
    formatted = format;
    # handle 0-360 vs. -180 to 180 for longitudes
    if not isLatitude:
        tval = value;
        if (use360):
            # taken from ucar.visad.GeoUtils.normalizeLongitude360
            # not sure why it's 361 instead of 360
            while ((tval < 0.) or (tval > 361.)):
                tval = 180. + math.remainder(tval - 180., 360.0);
        else:
            tval = normalize_longitude(tval);
        value = tval;
    
    pvalue = abs(value);
    # TODO:  Simplify this
    # convert to seconds then get the integer deg, min, seconds
    j        = int(3600.0 * pvalue);
    idegrees = int(j / 3600);
    decidegrees = (pvalue - int(pvalue));
    dminutes    = decidegrees * 60.;
    minutes     = int(dminutes);
    deciminutes = dminutes - minutes;
    dseconds    = deciminutes * 60.;
    seconds     = int(dseconds);
    deciseconds = dseconds - seconds;

    formatted = replace_decimal_portion(formatted, "d", decidegrees);
    formatted = replace_decimal_portion(formatted, "m", deciminutes);
    formatted = replace_decimal_portion(formatted, "s", deciseconds);
    formatted = formatted.replace("DD", str(idegrees));
    formatted = formatted.replace("MM",str(minutes).zfill(2));
    if 's' in format:
        formatted = formatted.replace("SS", str(seconds))
    else:
        formatted = formatted.replace("SS", str(int(round(dseconds))).zfill(2));

    #if (format.indexOf("H") >= 0):
    if 'H' in format:
        if (use360):                         # should we ignore or add E?
            formatted = formatted.replace("H", "",1);
        elif (abs(value) == 180):  # 180 line - subject to debate
            formatted = formatted.replace("H", "",1);
        elif (value < 0):  # South/West
            if isLatitude:
               formatted = formatted.replace("H","S",1)
            else:
               formatted = formatted.replace("H","W",1)
        elif (value > 0):  # North/East
            if isLatitude:
               formatted = formatted.replace("H","N",1)
            else:
               formatted = formatted.replace("H","E",1)
        else:                 # 0 line - subject to debate
            formatted = formatted.replace("H", "",1);
    elif ((value < 0) and (not use360)):
        formatted = "-" + formatted;
    
    return formatted.strip();

"""
 * Replace the decimal portion of a format string with the value
 *
 *  format  the format
 *  letter  the letter to replace (1 or more instances)
 *  decimalvalue  the value to replace it with
 *
 * return  the format with the appropriate value filled in
"""
def replace_decimal_portion(format, letter, decimalvalue):
    matches = re.findall(letter+'+',format)
    for match in matches:
        numchars = len(match)
        valueStr = str(round(decimalvalue * pow(10, numchars)));
        # account for zero;
        valueStr = valueStr.ljust(numchars, '0');
        #format   = matcher.replaceFirst(valueStr);
        format = format.replace(match,valueStr,1)
        break;
    
    return format;

"""
 * Normalize a longitude value to the range between -180 and 180.
 *
 * lonValue  longitude value to adjust (in degrees)
 * return adjusted value.
"""
def normalize_longitude(lonValue):
    while ((lonValue < -180.) or (lonValue > 180.)):
        if (lonValue < -180):
           lonValue = lonValue + 360.
        else:
           lonValue = lonValue - 360.;
    return lonValue;

