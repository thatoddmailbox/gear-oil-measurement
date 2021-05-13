import nanovnacmd as nv
import numpy as np
import matplotlib.pyplot as plt

from typing import List

from NanoVNASaver.Calibration import Calibration
from NanoVNASaver.RFTools import Datapoint, reflection_coefficient
from NanoVNASaver.Touchstone import Options, Touchstone

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

ts = Touchstone("./can3_oil_aluminum_fine.s1p")
ts.load()

f = np.array([d.freq for d in ts.s11data])
s11 = np.array([d.z for d in ts.s11data])

plt.plot(f/1e9,20*np.log10(np.abs(s11)),label='S11')
# plt.plot(f/1e9,20*np.log10(np.abs(s21)),label='S21')
plt.xlabel('Frequency (GHz)')
plt.ylabel('S parameter (dB)')
plt.ylim(-80,5)
plt.show()

plt.plot(f/1e9,np.angle(s11),label='S11')
# plt.plot(f/1e9,np.angle(s21),label='S21')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Phase (rad)')
plt.ylim(-np.pi,np.pi)
# plt.show()

resonant_frequency = calculate_resonance_frequency(ts)
print("resonant_frequency", resonant_frequency)

resonant_s11 = closest_point(ts.s11data, resonant_frequency)
resonant_gain = resonant_s11.gain
print("resonant gain:", resonant_gain)

half_power = resonant_gain + 3
print("half power:", half_power)
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
print(resonant_frequency)
print(intercepts)

first_intercept = None
second_intercept = None
for intercept in intercepts:
	if intercept > resonant_frequency:
		second_intercept = intercept
		break
	first_intercept = intercept
half_power_bandwidth = second_intercept - first_intercept
print(first_intercept, second_intercept, "half power bandwidth:", half_power_bandwidth)

# Q_loaded = half power bandwidth
q_loaded = resonant_frequency / half_power_bandwidth
print("q_loaded:", q_loaded)

resonance_band_from = 2.500 * 1e9
resonance_band_to = 2.870 * 1e9
# https://media.proquest.com/media/hms/ORIG/1/7ZvgC?hl=&cit%3Aauth=Shahid%2C+S%3BBall%2C+J+A+R%3BWells%2C+C+G%3BWen%2C+P&cit%3Atitle=Reflection+type+Q-factor+measurement+using+standard+least+squares+methods&cit%3Apub=IEE+proceedings.+Microwaves%2C+antennas+and+propagation.&cit%3Avol=5&cit%3Aiss=4&cit%3Apg=426&cit%3Adate=Mar+2011&ic=true&cit%3Aprod=Advanced+Technologies+%26+Aerospace+Collection&_a=ChgyMDIxMDUxMDIwMjIwMjQzNjo2Nzk2MDESBjEwMDg2NRoKT05FX1NFQVJDSCIMMTguMTguMjM5LjE2KgcxOTM2MzYxMgoxNjM2NzAyNTA1Og1Eb2N1bWVudEltYWdlQgEwUgZPbmxpbmVaAkZUYgNQRlRqCjIwMTEvMDMvMDFyCjIwMTEvMDMvMzF6AIIBKVAtMTAwNzg1Mi0xMjQ5Mi1DVVNUT01FUi0xMDAwMDIzMy01MDQ4NjM5kgEGT25saW5lygFMTW96aWxsYS81LjAgKFgxMTsgVWJ1bnR1OyBMaW51eCB4ODZfNjQ7IHJ2Ojg3LjApIEdlY2tvLzIwMTAwMTAxIEZpcmVmb3gvODcuMNIBElNjaG9sYXJseSBKb3VybmFsc5oCB1ByZVBhaWSqAi5PUzpFTVMtRG9jVmlld1BkZlVybFNlcnZpY2UtZ2V0TWVkaWFVcmxGb3JJdGVtygIPQXJ0aWNsZXxGZWF0dXJl0gIBWeICAU7yAgD6AgFOggMDV2ViigMcQ0lEOjIwMjEwNTEwMjAyMjAyNDMyOjg1Mzg4Mg%3D%3D&_s=3xxU4Nz07bp9dqMrWeemRkVM%2FxI%3D
# now that we have q, find gamma L, the closest to origin
gamma_L = None
gamma_L_dist = None
for p in ts.s11data:
	if p.freq < resonance_band_from or p.freq > resonance_band_to:
		continue
	gamma = p.z
	distance = abs(gamma)
	if not gamma_L_dist or distance < gamma_L_dist:
		# print(p, distance)
		gamma_L = gamma
		gamma_L_dist = distance
print("gamma_L:", gamma_L)

gamma_L_angle = np.angle(gamma_L, deg=False)
print("gamma_L angle (radians):", gamma_L_angle)
print("gamma_L angle (degrees):", np.rad2deg(gamma_L_angle))

diameter = 1 - abs(gamma_L)
print("diameter:", diameter)

kappa = diameter / (2 - diameter)
print("kappa:", kappa)

q_unloaded = q_loaded*(1 + kappa)
print("q_unloaded:", q_unloaded)