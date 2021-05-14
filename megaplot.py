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

for i in range(0, 8):
	for run in range(0, 2):
		ts = Touchstone("./data_water/water_"+str(i)+"_"+str(run)+".s1p")
		ts.load()

		f = np.array([d.freq for d in ts.s11data])
		s11 = np.array([d.z for d in ts.s11data])

		plt.plot(f/1e9,20*np.log10(np.abs(s11)),label='S11 ' + str(i) + "_" + str(run))
plt.xlabel('Frequency (GHz)')
plt.ylabel('S parameter (dB)')
plt.ylim(-80,5)
plt.legend()
plt.show()

# plt.plot(f/1e9,np.angle(s11),label='S11')
# # plt.plot(f/1e9,np.angle(s21),label='S21')
# plt.xlabel('Frequency (GHz)')
# plt.ylabel('Phase (rad)')
# plt.ylim(-np.pi,np.pi)
# # plt.show()