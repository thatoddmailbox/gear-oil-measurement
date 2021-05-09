### setup ###
import nanovnacmd as nv
import numpy as np
import matplotlib.pyplot as plt
#connect to the nanovna
vna = nv.connect()

### calibration ###
#calculate the frequency range for the tdr lab
#use a frequency range with 101 steps (102 when including the extrapolated DC value)
fwithDC = np.linspace(1e9,4e9,102)

#set the frequencies
vna.setSweep(fwithDC[1],fwithDC[-1])

# #2 port calibration - short, open, load, through, isolation (cap the two lines w/ 50 ohms)
calibration = nv.calibrate2port(vna)

calibration.save("calibration.txt")