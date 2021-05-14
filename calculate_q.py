import nanovnacmd as nv
import numpy as np
import matplotlib.pyplot as plt

from typing import List

from NanoVNASaver.Calibration import Calibration
from NanoVNASaver.RFTools import Datapoint, reflection_coefficient
from NanoVNASaver.Touchstone import Options, Touchstone

import sys
sys.path.append("../pySmithPlot/")
from smithplot import SmithAxes

from scipy import optimize
from numpy import sqrt
from math import pi, sin, cos

water_ranges = {
    './data_water/water_0_0.s1p': (2.518072e9, 3.03012e9),
    './data_water/water_0_1.s1p': (2.524096e9, 3.03012e9),
    './data_water/water_1_0.s1p': (2.572289e9, 3.138554e9),
    './data_water/water_2_0.s1p': (2.596386e9, 3.150602e9),
    './data_water/water_2_1.s1p': (2.603225806e9, 3.158064516e9),
    './data_water/water_3_0.s1p': (2.719354839e9, 3.156627e9),
    './data_water/water_3_1.s1p': (2.609677419e9, 2.987096774e9),
    './data_water/water_4_0.s1p': (2.706451613e9, 3.156627e9),
    './data_water/water_4_1.s1p': (2.680645161e9, 3.151612903e9),
    './data_water/water_5_0.s1p': (2.761290323e9, 3.258064516e9),
    './data_water/water_5_1.s1p': (2.759036e9, 3.264516129e9),
    './data_water/water_6_0.s1p': (2.809677419e9, 3.329032258e9),
    './data_water/water_6_1.s1p': (2.835483871e9, 3.331325e9),
    './data_water/water_7_0.s1p': (2.861290323e9, 3.331325e9),
    './data_water/water_7_1.s1p': (2.883870968e9, 3.358064516e9),
}

def norm(x: float, y: float) -> List[float]:
	z = x + 1j*y
	result = (1 + z)/(1 - z)
	return result.real, result.imag

def plot_s11(ts: Touchstone, circle: (float, float, float), trim_start: List[Datapoint], trim_end: List[Datapoint]):
	f = np.array([d.freq for d in ts.s11data])
	s11 = np.array([d.z for d in ts.s11data])

	trim_data = trim_start+trim_end
	trim_f = np.array([d.freq for d in trim_data])
	trim_s11 = np.array([d.z for d in trim_data])

	plt.plot(f/1e9,20*np.log10(np.abs(s11)),label='S11')
	# plt.plot(f/1e9,20*np.log10(np.abs(s21)),label='S21')
	plt.xlabel('Frequency (GHz)')
	plt.ylabel('S parameter (dB)')
	plt.ylim(-80,5)
	# plt.show()

	# plt.plot(f/1e9,np.angle(s11),label='S11')
	# # plt.plot(f/1e9,np.angle(s21),label='S21')
	# plt.xlabel('Frequency (GHz)')
	# plt.ylabel('Phase (rad)')
	# plt.ylim(-np.pi,np.pi)
	# # plt.show()

	plt.figure(figsize=(6, 6))

	xc_2, yc_2, R_2 = circle
	print(xc_2, yc_2)
	print(R_2)
	circle_points = []
	for i in range(0, 360):
		x = R_2 * cos(i * (pi / 180))
		y = R_2 * sin(i * (pi / 180))
		x += xc_2
		y += yc_2
		circle_points.append(x + y*1j)

	ax = plt.subplot(1, 1, 1, projection='smith')
	# plt.plot([10, 100], markevery=1)

	plt.plot(s11, label="default", datatype=SmithAxes.S_PARAMETER)
	plt.plot(trim_s11, label="default", datatype=SmithAxes.S_PARAMETER)
	# plt.plot(xc_2 + 1j*yc_2)
	plt.plot(circle_points, linewidth=0.25)

	# leg = plt.legend(loc="lower right", fontsize=12)
	plt.title("Smith Chart")

	plt.show()

def closest_point(data: List[Datapoint], freq: int) -> Datapoint:
	closest = None
	closest_dist = None
	for point in data:
		if not closest_dist or abs(point.freq - freq) < closest_dist:
			closest = point
			closest_dist = abs(point.freq - freq)
	return closest

def calculate_resonance_frequency(ts: Touchstone) -> float:
	lowest_s11_freq = None
	lowest_s11 = None
	for p in ts.s11data:
		# print(p.freq, p.gain)
		if not lowest_s11 or abs(p.gain) > lowest_s11:
			lowest_s11_freq = p.freq
			lowest_s11 = abs(p.gain)
	return lowest_s11_freq

def calculate_half_power_bandwidth(ts: Touchstone, resonant_frequency: float, resonant_gain: float) -> float:
	# find half power value in decibels
	half_power = resonant_gain + 3
	print("half power:", half_power)

	# find intercepts at half power value
	first = True
	above = False
	intercepts = []
	for p in ts.s11data:
		value = 20*np.log10(np.abs(p.z))
		if first:
			first = False
			above = value > half_power
			continue
		if (above and value < half_power) or (not above and value > half_power):
			above = not above
			intercepts.append(p.freq)

	assert(len(intercepts) >= 2)
	print("resonant_frequency", resonant_frequency)
	print(intercepts)

	# find intercepts of half power
	first_intercept = None
	second_intercept = None
	for intercept in intercepts:
		if intercept > resonant_frequency:
			second_intercept = intercept
			break
		first_intercept = intercept
	print(first_intercept, second_intercept)
	half_power_bandwidth = second_intercept - first_intercept
	print(first_intercept, second_intercept, "half power bandwidth:", half_power_bandwidth)

	return half_power_bandwidth

def fit_circle(ts: Touchstone) -> (float, float, float):
	s11_re = np.array([d.re for d in ts.s11data])
	s11_im = np.array([d.im for d in ts.s11data])
	# print(s11_re)

	def calc_R(xc, yc):
		""" calculate the distance of each 2D points from the center (xc, yc) """
		return sqrt((s11_re-xc)**2 + (s11_im-yc)**2)

	def f_2(c):
		""" calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
		Ri = calc_R(*c)
		return Ri - Ri.mean()

	center_estimate = 0, 0
	center_2, ier = optimize.leastsq(f_2, center_estimate)

	xc_2, yc_2 = center_2
	Ri_2       = calc_R(*center_2)
	R_2        = Ri_2.mean()
	residu_2   = sum((Ri_2 - R_2)**2)

	return (xc_2, yc_2, R_2)

filename = "./data_water/water_7_1.s1p"
ts = Touchstone(filename)
ts.load()

resonance_band_from = 2.500 * 1e9
resonance_band_to = 3.0 * 1e9
if filename in water_ranges:
	resonance_band_from, resonance_band_to = water_ranges[filename]
print(resonance_band_from, resonance_band_to)

found_from = False
index_from = None
index_to = None
for i, p in enumerate(ts.s11data):
	if not found_from:
		if p.freq < resonance_band_from:
			continue
		index_from = i
		found_from = True
	else:
		if p.freq < resonance_band_to:
			continue
		index_to = i
		break
trim_start = ts.s11data[:index_from+1]
trim_end = ts.s11data[index_to:]
ts.s11data = ts.s11data[index_from:index_to+1]

circle = fit_circle(ts)
plot_s11(ts, circle, trim_start, trim_end)

resonant_frequency = calculate_resonance_frequency(ts)
print("resonant_frequency", resonant_frequency)

resonant_s11 = closest_point(ts.s11data, resonant_frequency)
resonant_gain = resonant_s11.gain
print("resonant gain:", resonant_gain)

half_power_bandwidth = calculate_half_power_bandwidth(ts, resonant_frequency, resonant_gain)

# Q_loaded = half power bandwidth
q_loaded = resonant_frequency / half_power_bandwidth
print("q_loaded:", q_loaded)

# resonance_band_from = 2.500 * 1e9
# resonance_band_to = 2.870 * 1e9
# https://media.proquest.com/media/hms/ORIG/1/7ZvgC?hl=&cit%3Aauth=Shahid%2C+S%3BBall%2C+J+A+R%3BWells%2C+C+G%3BWen%2C+P&cit%3Atitle=Reflection+type+Q-factor+measurement+using+standard+least+squares+methods&cit%3Apub=IEE+proceedings.+Microwaves%2C+antennas+and+propagation.&cit%3Avol=5&cit%3Aiss=4&cit%3Apg=426&cit%3Adate=Mar+2011&ic=true&cit%3Aprod=Advanced+Technologies+%26+Aerospace+Collection&_a=ChgyMDIxMDUxMDIwMjIwMjQzNjo2Nzk2MDESBjEwMDg2NRoKT05FX1NFQVJDSCIMMTguMTguMjM5LjE2KgcxOTM2MzYxMgoxNjM2NzAyNTA1Og1Eb2N1bWVudEltYWdlQgEwUgZPbmxpbmVaAkZUYgNQRlRqCjIwMTEvMDMvMDFyCjIwMTEvMDMvMzF6AIIBKVAtMTAwNzg1Mi0xMjQ5Mi1DVVNUT01FUi0xMDAwMDIzMy01MDQ4NjM5kgEGT25saW5lygFMTW96aWxsYS81LjAgKFgxMTsgVWJ1bnR1OyBMaW51eCB4ODZfNjQ7IHJ2Ojg3LjApIEdlY2tvLzIwMTAwMTAxIEZpcmVmb3gvODcuMNIBElNjaG9sYXJseSBKb3VybmFsc5oCB1ByZVBhaWSqAi5PUzpFTVMtRG9jVmlld1BkZlVybFNlcnZpY2UtZ2V0TWVkaWFVcmxGb3JJdGVtygIPQXJ0aWNsZXxGZWF0dXJl0gIBWeICAU7yAgD6AgFOggMDV2ViigMcQ0lEOjIwMjEwNTEwMjAyMjAyNDMyOjg1Mzg4Mg%3D%3D&_s=3xxU4Nz07bp9dqMrWeemRkVM%2FxI%3D
# now that we have q, find gamma L, the closest to origin
# gamma_L = None
# gamma_L_dist = None
# for p in ts.s11data:
# 	if p.freq < resonance_band_from or p.freq > resonance_band_to:
# 		continue
# 	gamma = p.z
# 	distance = abs(gamma)
# 	if not gamma_L_dist or distance < gamma_L_dist:
# 		# print(p, distance)
# 		gamma_L = gamma
# 		gamma_L_dist = distance
# print("gamma_L:", gamma_L)

# gamma_L_angle = np.angle(gamma_L, deg=False)
# print("gamma_L angle (radians):", gamma_L_angle)
# print("gamma_L angle (degrees):", np.rad2deg(gamma_L_angle))

diameter = 2*circle[2]
# diameter = 1 - abs(gamma_L)
print("diameter:", diameter)

kappa = diameter / (2 - diameter)
print("kappa:", kappa)

q_unloaded = q_loaded*(1 + kappa)
print("q_unloaded:", q_unloaded)