### setup ###
import nanovnacmd as nv
import numpy as np
import matplotlib.pyplot as plt

from NanoVNASaver.Calibration import Calibration

# load calibration
calibration = Calibration()
calibration.load("calibration.txt")

# connect to nanovna
vna = nv.connect()

fwithDC = np.linspace(0,1.01e9,102)

vna.setSweep(fwithDC[1],fwithDC[-1])

f,s11,s21 = nv.measure(vna, calibration)

plt.plot(f/1e9,20*np.log10(np.abs(s11)),label='S11')
plt.plot(f/1e9,20*np.log10(np.abs(s21)),label='S21')
plt.xlabel('Frequency (GHz)')
plt.ylabel('S parameter (dB)')
plt.ylim(-80,5)
plt.show()

plt.plot(f/1e9,np.angle(s11),label='S11')
plt.plot(f/1e9,np.angle(s21),label='S21')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Phase (rad)')
plt.ylim(-np.pi,np.pi)
plt.show()