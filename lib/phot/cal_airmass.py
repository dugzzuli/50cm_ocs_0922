import numpy as np
from scipy import integrate
from astropy.time import Time
from astropy.time import TimeDelta
from astropy import coordinates as coord, units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import os, warnings, sys
from loguru import logger

###defining function###
def airmass(hour_angle,latitude,dec):
	"""
	airmass at hour angle t
	"""
	x = 1/(np.sin(latitude)*np.sin(dec)+np.cos(latitude)*np.cos(dec)*np.cos(hour_angle))  #t is hour angle
	delta_x = 0.00186*(x-1)+0.002875*((x-1)**2)+0.0008083*((x-1)**3)
	airmass = x - delta_x
	return float(airmass)

def cam_old(obstime,exp_time,ra,dec):
	ra = float(ra)*(np.pi/180)
	dec = float(dec)*(np.pi/180)

	###calculating the sidereal time at the beginning of observation###
	lon_deg,lon_min,lon_sec = 100,1,48
	longitude = (float(lon_deg)+float(lon_min)/60+float(lon_sec)/3600)
	lat_d,lat_arcm,lat_arcs = 26,42,32
	latitude = (float(lat_d)+float(lat_arcm)/60+float(lat_arcs)/3600)*(np.pi/180)
	time_standard = 'utc'
	
	#the jd/mgst/lmst/tdb at the beginning of exposure
	#obstime = '{}-{}-{}T{}:{}:{}'.format(int(year), int(month), int(day), int(hour), int(minute), float(second))
	scale = time_standard
	location = ('{}d'.format(str(longitude)),'{}d'.format(str(latitude/(np.pi/180))))
	obstime_start = Time(obstime,format="isot",scale=scale,location=location)
	
	ut1_start = obstime_start.ut1.iso
	jd_start = obstime_start.ut1.jd
	gmst_start = obstime_start.ut1.sidereal_time('mean', 'greenwich','IAU2006')
	lmst_start = obstime_start.ut1.sidereal_time('mean',model='IAU2006')
	tdb_start = obstime_start.ut1.tdb.iso
	
	#the jd/mgst/lmst/tdb at the mid-point of exposure
	mean_second = exp_time/2
	delta_time_mid = TimeDelta(mean_second,format='sec')
	obstime_mid = obstime_start + delta_time_mid
	
	ut1_mid = obstime_mid.ut1.iso
	jd_mid = obstime_mid.ut1.jd
	gmst_mid = obstime_mid.ut1.sidereal_time('mean', 'greenwich','IAU2006')
	lmst_mid = obstime_mid.ut1.sidereal_time('mean',model='IAU2006')
	tdb_mid = obstime_mid.ut1.tdb.iso
	
	###the hour angle (equal time interval)###
	n = 5
	lmst_hour = np.array(lmst_start)*15*(np.pi/180)
	t0 = lmst_hour - ra  # initial hour angle h
	hour_angle = []
	for i in np.arange(0,(float(exp_time)+float(exp_time)/n),(float(exp_time)/n)):
		t = t0 + (i*15/3600)*(np.pi/180)
		hour_angle.append(t)
	
	###method 1:calculating the airmass by chosing the average time point###
	t_start = hour_angle[0]     #the beginning point of hour angle of observation
	t_end = hour_angle[len(hour_angle)-1]    #the ending point of hour angle of observation
	ave_t1 = (t_end+t_start)/2
	
	airmass_start = airmass(t_start,latitude,dec)
	airmass_end = airmass(t_end,latitude,dec)
	airmass_mean = airmass(ave_t1,latitude,dec)
	ave_airmass = round((airmass_start + 4*airmass_mean + airmass_end)/6,4)
	return(ave_airmass,jd_mid,tdb_mid)

def cam(obstime,exp_time,ra,dec):
	# print(obstime)
	# print(exp_time)
	ra = float(ra)*(np.pi/180)
	dec = float(dec)*(np.pi/180)
	###calculating the sidereal time at the beginning of observation###
	lon_deg,lon_min,lon_sec = 100,1,48
	longitude = (float(lon_deg)+float(lon_min)/60+float(lon_sec)/3600)
	lat_d,lat_arcm,lat_arcs = 26,42,32
	latitude = (float(lat_d)+float(lat_arcm)/60+float(lat_arcs)/3600)*(np.pi/180)
	time_standard = 'utc'
	#the jd/mgst/lmst/tdb at the beginning of exposure
	#obstime = '{}-{}-{}T{}:{}:{}'.format(int(year), int(month), int(day), int(hour), int(minute), float(second))
	scale = time_standard
	location = ('{}d'.format(str(longitude)),'{}d'.format(str(latitude/(np.pi/180))))
	obstime_start = Time(obstime,format="isot",scale=scale,location=location)
	#ut1_start = obstime_start.ut1.iso
	#jd_start = obstime_start.ut1.jd
	#gmst_start = obstime_start.ut1.sidereal_time('mean', 'greenwich','IAU2006')
	lmst_start = obstime_start.ut1.sidereal_time('mean',model='IAU2006')
	#tdb_start = obstime_start.ut1.tdb.iso
	#the jd/mgst/lmst/tdb at the mid-point of exposure
	mean_second = exp_time/2
	delta_time_mid = TimeDelta(mean_second,format='sec')
	obstime_mid = obstime_start + delta_time_mid
	#ut1_mid = obstime_mid.ut1.iso
	jd_mid = obstime_mid.ut1.jd
	#gmst_mid = obstime_mid.ut1.sidereal_time('mean', 'greenwich','IAU2006')
	#lmst_mid = obstime_mid.ut1.sidereal_time('mean',model='IAU2006')
	tdb_mid = obstime_mid.ut1.tdb.iso
	###the hour angle (equal time interval)###
	n = 5
	lmst_hour = np.array(lmst_start)*15*(np.pi/180)
	t0 = lmst_hour - ra  # initial hour angle h
	hour_angle = []
	for i in np.arange(0,(float(exp_time)+float(exp_time)/n),(float(exp_time)/n)):
		t = t0 + (i*15/3600)*(np.pi/180)
		hour_angle.append(t)
	###method 1:calculating the airmass by chosing the average time point###
	t_start = hour_angle[0]     #the beginning point of hour angle of observation
	t_end = hour_angle[len(hour_angle)-1]    #the ending point of hour angle of observation
	ave_t1 = (t_end+t_start)/2
	airmass_start = airmass(t_start,latitude,dec)
	airmass_end = airmass(t_end,latitude,dec)
	airmass_mean = airmass(ave_t1,latitude,dec)
	ave_airmass = round((airmass_start + 4*airmass_mean + airmass_end)/6,4)
	return(ave_airmass,jd_mid,tdb_mid)
