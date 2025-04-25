import numpy as np

variables = { 
    'sfc_met' : [
        {   
            'varname' : 'pressure',
            'valid_range' : np.array([850., 1100.], np.float32),
            'attrs' : {
                'long_name' : 'Surface Air Pressure',
                'units' : 'hPa',
                'standard_name' : 'surface_air_pressure',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'temperature',
            'valid_range' : np.array([-90., 50.], np.float32),
            'attrs' : {
                'long_name' : 'Air Temperature',
                'units' : 'degC',
                'standard_name' : 'air_temperature',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'relative_humidity',
            'attrs' : {
                'long_name' : 'Relative Humidity',
                'units' : '%',
                'standard_name' : 'relative_humidity',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'wind_speed',
            'attrs' : {
                'long_name' : 'Scalar Wind Speed',
                'units' : 'm/s',
                'standard_name' : 'wind_speed',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'vecwind_speed',
            'attrs' : {
                'long_name' : 'Vector Wind Speed',
                'units' : 'm/s',
                'standard_name' : 'wind_speed',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'wind_direction',
            'attrs' : {
                'long_name' : 'Wind Direction',
                'units' : 'degrees',
                'standard_name' : 'wind_from_direction',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'stdwind_direction',
            'attrs' : {
                'long_name' : 'Wind Direction Standard Deviation',
                'units' : 'degrees',
                'standard_name' : 'wind_from_direction',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'maxwind_speed',
            'attrs' : {
                'long_name' : 'Maximum Wind Speed',
                'units' : 'm/s',
                'standard_name' : 'wind_speed_of_gust',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'battery_voltage',
            'attrs' : {
                'long_name' : 'Battery Voltage',
                'units' : 'volts',
                'standard_name' : 'battery_voltage',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'solar_radiation',
            'attrs' : {
                'long_name' : 'Solar Radiation',
                'units' : 'W m-2',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'net_radiation',
            'attrs' : {
                'long_name' : 'Net Radiation',
                'units' : 'W m-2',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'precipitation',
            'attrs' : {
                'long_name' : 'Precipitation',
                'units' : 'mm',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'soil_temperature_10cm',
            'attrs' : {
                'long_name' : 'Soil Temperature at 10cm',
                'units' : 'degC',
                'standard_name' : 'soil_temperature',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'soil_temperature_15cm',
            'attrs' : {
                'long_name' : 'Soil Temperature at 15cm',
                'units' : 'degC',
                'standard_name' : 'soil_temperature',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'soil_reflectometer_period_10cm',
            'attrs' : {
                'long_name' : 'Soil Reflectometer Output Period at 10cm',
                'units' : 'microsecond',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'soil_reflectometer_period_15cm',
            'attrs' : {
                'long_name' : 'Soil Reflectometer Output Period at 15cm',
                'units' : 'microsecond',
                '_FillValue' : 99999.,
                'missing_value' : 99999.
            }
        },
        {   
            'varname' : 'base_time',
            'attrs' : {
                'long_name' : 'Base Time',
                'units' : 'seconds since 1970-01-01 00:00:00'
            }
        },
        {   
            'varname' : 'time_offset',
            'attrs' : {
                'long_name' : 'Time Offset from Base Time',
                'units' : 'seconds'
            }
        },
        {   
            'varname' : 'time',
            'attrs' : {
                'long_name' : 'time',
                'units' : 'seconds since 1970-01-01 00:00:00',
            }
        },
    ],
    'sounding' : [
        {   
            'varname' : 'p1',
            'attrs' : {
                'long_name' : 'Air Pressure',
                'units' : 'hPa',
                'standard_name' : 'air_pressure',
            }
        },
        {   
            'varname' : 't',
            'attrs' : {
                'long_name' : 'Air Temperature',
                'units' : 'degC',
                'standard_name' : 'air_temperature',
            }
        },
        {   
            'varname' : 'q',
            'attrs' : {
                'long_name' : 'Specific Humidity',
                'units' : 'g/kg',
                'standard_name' : 'specific_humidity',
            }
        },
        {   
            'varname' : 'u',
            'attrs' : {
                'long_name' : 'Eastward Wind',
                'units' : 'm/s',
                'standard_name' : 'eastward_wind',
            }
        },
        {   
            'varname' : 'v',
            'attrs' : {
                'long_name' : 'Northward Wind',
                'units' : 'm/s',
                'standard_name' : 'northward_wind',
            }
        },
        {   
            'varname' : 'g',
            'attrs' : {
                'long_name' : 'Geopotential Height',
                'units' : 'm',
                'standard_name' : 'geopotential_height',
            }
        },
        {   
            'varname' : 'time',
            'attrs' : {
                'long_name' : 'time',
                'units' : 'seconds since 1970-01-01 00:00:00',
            }
        },
    ],
#    'height','wind_speed','wind_direction','qc','rv1','rv2','rv3','snr1','snr2','snr3','qc1','qc2','qc3'
    'rwp' : [
        {   
            'varname' : 'base_time',
            'attrs' : {
                'long_name' : 'Base Time',
                'units' : 'seconds since 1970-01-01 00:00:00'
            }
        },
        {   
            'varname' : 'time_offset',
            'attrs' : {
                'long_name' : 'Time Offset from Base Time',
                'units' : 'seconds'
            }
        },
        {   
            'varname' : 'time',
            'attrs' : {
                'long_name' : 'time',
                'units' : 'seconds since 1970-01-01 00:00:00',
            }
        },
        {   
            'varname' : 'height',
            'attrs' : { 'long_name' : 'Height above ground level',
                'units' : 'km',
                'valid_min' : 0.,
                'valid_max' : 8.,
                'missing_value' : 999999.,
                'standard_name' : "height",
                '_FillValue' : 999999.
             }
        },
        {   
            'varname' : 'wind_speed',
            'attrs' : {
                'long_name' : 'Wind Speed',
                'units' : 'm/s',
                'standard_name' : 'wind_speed',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'wind_direction',
            'attrs' : {
                'long_name' : 'Wind Direction',
                'units' : 'degrees',
                'standard_name' : 'wind_from_direction',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'wind_quality',
            'attrs' : {
                'long_name' : 'Quality Control For Wind',
                'flag_values'  : np.array([0, 2, 7, 8, 9]),
                'flag_meanings'  : "valid estimated suspect invalid missing",
                'standard_name' : 'quality_flag',
            }
        },
        {   
            'varname' : 'radial_velocity',
            'attrs' : {
                'long_name' : 'Radial Velocity',
                'units' : 'm s-1',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'radial_velocity_snr',
            'attrs' : {
                'long_name' : 'Signal To Noise Ratio',
                'units' : 'dB',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'radial_velocity_navg',
            'attrs' : {
                'long_name' : 'Number of Records in Radial Velocity Average',
                'units' : 1,
                '_FillValue' : 0,
                'missing_value' : 0
            }
        },
        {   
            'varname' : 'radial_velocity_quality',
            'attrs' : {
                'long_name' : 'Radial Velocity Average Quality',
                'units' : 1,
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'rv1',
            'attrs' : {
                'long_name' : 'Radial Velocity Beam 1',
                'units' : 'm s-1',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'rv2',
            'attrs' : {
                'long_name' : 'Radial Velocity Beam 2',
                'units' : 'm s-1',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'rv3',
            'attrs' : {
                'long_name' : 'Radial Velocity Beam 3',
                'units' : 'm s-1',
                '_FillValue' : 999999.,
                'missing_value' : 999999.
            }
        },
        {   
            'varname' : 'snr1',
            'attrs' : {
                'long_name' : 'Signal To Noise Ratio Beam 1',
                'units' : 'dB',
                '_FillValue' : 999999.,
                'missing_value' : 999999.,
            }
        },
        {   
            'varname' : 'snr2',
            'attrs' : {
                'long_name' : 'Signal To Noise Ratio Beam 2',
                'units' : 'dB',
                '_FillValue' : 999999.,
                'missing_value' : 999999.,
            }
        },
        {   
            'varname' : 'snr3',
            'attrs' : {
                'long_name' : 'Signal To Noise Ratio Beam 3',
                'units' : 'dB',
                '_FillValue' : 999999.,
                'missing_value' : 999999.,
            }
        },
        {   
            'varname' : 'qc1',
            'attrs' : {
                'long_name' : 'Quality Control Radial Velocity Beam 1',
            }
        },
        {   
            'varname' : 'qc2',
            'attrs' : {
                'long_name' : 'Quality Control Radial Velocity Beam 2',
            }
        },
        {   
            'varname' : 'qc3',
            'attrs' : {
                'long_name' : 'Quality Control Radial Velocity Beam 3',
            }
        },
        {   
            'varname' : 'wind_u',
            'attrs' : {
                'long_name' : 'Eastward Wind',
                'units' : 'm/s',
                'standard_name' : 'eastward_wind',
                '_FillValue' : 999999.,
                'missing_value' : 999999.,
            }
        },
        {   
            'varname' : 'wind_v',
            'attrs' : {
                'long_name' : 'Northward Wind',
                'units' : 'm/s',
                'standard_name' : 'northward_wind',
                '_FillValue' : 999999.,
                'missing_value' : 999999.,
            }
        },
#    'height','tv_uc','tv','w','qc_tv_uc','qc_tv','qc_w','cnt_tv_uc','cnt_tv','cnt_w','snr_tv_uc','snr_tv','snr_w']
        {   
            'varname' : 'tv_uc',
            'attrs' : {
                'long_name' : 'RASS Temperature (uncorrected)',
                'units' : 'degC',
                'standard_name' : 'air_temperature',
            }
        },
        {   
            'varname' : 'tv_c',
            'attrs' : {
                'long_name' : 'RASS Temperature (corrected)',
                'units' : 'degC',
                'standard_name' : 'air_temperature',
            }
        },
        {   
            'varname' : 'w',
            'attrs' : {
                'long_name' : 'Vertical Wind',
                'units' : 'm s-1',
                'standard_name' : 'upward_air_velocity',
            }
        },
        {   
            'varname' : 'qc_tv_uc',
            'attrs' : {
                'long_name' : 'Quality Control Flag for RASS Temperature (uncorrected)',
                'missing_value' : 9.,
                'flag_values'  : np.array([0, 9]),
                'flag_meanings'  : "valid missing",
                'standard_name' : 'quality_flag',
            }
        },
        {   
            'varname' : 'qc_tv_c',
            'attrs' : {
                'long_name' : 'Quality Control Flag for RASS Temperature (corrected)',
                'missing_value' : 9.,
                'flag_values'  : np.array([0, 9]),
                'flag_meanings'  : "valid missing",
                'standard_name' : 'quality_flag',
            }
        },
        {   
            'varname' : 'qc_w',
            'attrs' : {
                'long_name' : 'Quality Control Flag for Vertical Wind',
                'missing_value' : 9.,
                'flag_values'  : np.array([0, 9]),
                'flag_meanings'  : "valid missing",
                'standard_name' : 'quality_flag',
            }
        },
        {   
            'varname' : 'cnt_tv_uc',
            'attrs' : {
                'long_name' : 'Number of Records for Average of RASS Temperature (uncorrected)',
            }
        },
        {   
            'varname' : 'cnt_tv_c',
            'attrs' : {
                'long_name' : 'Number of Records for Average of RASS Temperature (corrected)',
            }
        },
        {   
            'varname' : 'cnt_w',
            'attrs' : {
                'long_name' : 'Number of Records for Average of Vertical Wind',
            }
        },
        {   
            'varname' : 'snr_tv_uc',
            'attrs' : {
                'long_name' : 'Average Signal-to-Noise (SNR) of RASS Temperature (uncorrected)',
                'units' : 'dB',
            }
        },
        {   
            'varname' : 'snr_tv_c',
            'attrs' : {
                'long_name' : 'Average Signal-to-Noise (SNR) of RASS Temperature (corrected)',
                'units' : 'dB',
            }
        },
        {   
            'varname' : 'snr_w',
            'attrs' : {
                'long_name' : 'Average Signal-to-Noise (SNR) of Vertical Wind',
                'units' : 'dB',
            }
        },
    ],
     'disdrometer' : [
         {
             'varname' : 'base_time',
             'attrs' : {
                 'long_name' : 'Base Time',
                 'units' : 'seconds since 1970-01-01 00:00:00'
             }
         },
         {
             'varname' : 'time_offset',
             'attrs' : {
                 'long_name' : 'Time Offset from Base Time',
                 'units' : 'seconds'
             }
         },
         {
             'varname' : 'time_interval',
             'attrs' : {
                 'long_name' : 'Observation Time Interval',
                 'units' : 'seconds'
             }
         },
         {
             'varname' : 'time',
             'attrs' : {
                 'long_name' : 'time',
                 'units' : 'seconds since 1970-01-01 00:00:00',
              }
         },
         {
             'varname' : 'Rate',
             'attrs' : {
                 'long_name' : 'Precipitation Rate',
                 'units' : 'mm/h',
              }
         },
         {
             'varname' : 'Amount',
             'attrs' : {
                 'long_name' : 'Precipitation Amount',
                 'units' : 'mm',
              }
         },
         {
             'varname' : 'AmountSum',
             'attrs' : {
                 'long_name' : 'Event Precipitation Amount',
                 'units' : 'mm',
              }
         },
    ],
}
