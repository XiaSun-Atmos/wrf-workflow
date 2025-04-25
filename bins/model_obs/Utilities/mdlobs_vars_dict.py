import numpy as np

variables = { 
    'mdl_met' : [
        {   
            'varname' : 'precipitation',
            'attrs' : {
                'long_name' : 'Precipitation',
                'standard_name' : 'precipitation_amount',
                'units' : 'mm',
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
            'varname' : 'time',
            'attrs' : {
                'long_name' : 'time',
                'units' : 'seconds since 1970-01-01 00:00:00',
            }
        },
        {   
            'varname' : 'forecast',
            'attrs' : {
                'long_name' : 'forecast hour',
                'units' : 'hour',
            }
        },
    ],

}
