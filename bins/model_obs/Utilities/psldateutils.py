from datetime import datetime,timezone,timedelta


# String can be yyjjj or yyjjjhh
def psldatestr_to_datetime(psldate):
    hh=0
    year=int(psldate[0:2])
    jday=int(psldate[2:5])
    if len(psldate) > 5:
      hh=int(psldate[5:])
    return psldate_to_datetime(year,jday,hour=hh)

    

def psldate_to_datetime(year,jday,hour=0,min=0):
  yr=int(year)
  if (yr<100):
      yr=2000+yr
  jd=int(jday)
  hrstr=str(hour)
  if len(hrstr) == 4:
    hr=int(hrstr[0:2])
    mn=int(hrstr[2:])
  else:
    hr=int(hrstr)
    mn=int(min)
  dt=datetime(yr,1,1,tzinfo=timezone.utc)+timedelta(days=(jd-1),minutes=mn,hours=hr)
  return dt

# return the datetime as a year,jday tuple
def datetime_to_yyyyddd(dt):
  tt=dt.timetuple();
  return [tt.tm_year,tt.tm_yday]


# datestring as yyyymmdd (e.g., 20230301) to yyyyjjj  (2023060)
def ymd_to_yj(datestr):
   dt=datetime(int(datestr[0:4]),int(datestr[4:6]),int(datestr[6:8]),tzinfo=timezone.utc)
   tt=datetime_to_yyyyddd(dt)
   return str(tt[0])+str(tt[1]).zfill(3)
  
