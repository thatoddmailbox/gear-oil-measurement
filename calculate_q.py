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

ts = Touchstone("./test1_good.s1p")
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

resonant_frequency = 2.879 * 1e9
resonant_s11 = closest_point(ts.s11data, resonant_frequency)

half_power = -3 # dB
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
print(q_loaded)

# now that we have q, find gamma L, the closest to origin
gamma_L = None
gamma_L_dist = None
for p in ts.s11data:
	gamma = reflection_coefficient(p.z)
	distance = abs(gamma)
	if not gamma_L_dist or distance > gamma_L_dist:
		gamma_L = gamma
		gamma_L_dist = distance
print(gamma_L)