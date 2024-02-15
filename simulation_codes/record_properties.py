import numpy as np
from sunRay import plasmaFreq as pfreq
from sunRay import densityModel as dm
from sunRay import scattering as scat 
#from sunRay import showPlot as SP   # matplotlib will slow down the run
from sunRay.parameters import c,c_r,R_S  # physics parameters
from sunRay.parameters import dev_u  # computation device
import sunRay.statisticalRays as raystat
import torch
import time
from tqdm import tqdm # for processing bar
import datetime



def record_properties2(dist, dict_dist, rr_cur, rx_cur, ry_cur, rz_cur, t_current):
	n = torch.where(torch.isnan(rr_cur) == False)[0]   #sometimes there are nan values that appear  (Note this line slows down the code by amout 10 it/s)
	#print(dist, torch.max(rr_cur[n]), torch.min(rr_cur[n]))
	if (torch.max(rr_cur[n]) >= dist) and (torch.min(rr_cur[n]) < dist):
		#start = time.process_time()
		q = torch.where((rr_cur >= dist) & (dict_dist['time'] == 0))[0]
		if q.size()[0] > 0:
			dict_dist['time'][q] = t_current.cpu()
			dict_dist['x'][q] = rx_cur[q]
			dict_dist['y'][q] = ry_cur[q]
			dict_dist['z'][q] = rz_cur[q]

		#print('Process:',time.process_time() - start)

	return dict_dist










def record_properties_thetarng(dist, t_dist_theta, theta_low, theta_up, rr_cur, rx_cur, ry_cur, rz_cur, t_current):
	#angles are in degrees
	R, x, y, z = rr_cur.cpu().numpy(), rx_cur.cpu().numpy(),  ry_cur.cpu().numpy(),  rz_cur.cpu().numpy()
	

	theta = np.arccos(z / np.sqrt(x**2 + z**2)) * 180 / np.pi


	i = np.where((R > dist*0.999) & (R < dist*1.001) & (theta > theta_low) & (theta < theta_up))[0]    #which photons are very close to one 1AU (i.e dist)
	if len(i) > 0:
		j = np.where(t_dist_theta != 0)[0]						 #which photons have already passed once through 1AU
		ij = np.concatenate((i, j), axis = 0)
		u, c = np.unique(ij, return_counts = True)
		dupl = u[c > 1]  			  #these are the photons which have already passed through 1AU (i.e duplicates)
		q = np.setdiff1d(i, dupl)     #gets rid of the photons  in array iwhich have already passed through 1AU

		t_dist_theta[q] = t_current

	return t_dist_theta