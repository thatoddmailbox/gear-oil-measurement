### setup ###
import nanovnacmd as nv
import numpy as np
import matplotlib.pyplot as plt

from NanoVNASaver.Calibration import Calibration
from NanoVNASaver.Touchstone import Options, Touchstone

ts = Touchstone("./test1_good.s1p")
ts.load()

print(ts.opts)

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
plt.show()

resonant_frequency = 2.879