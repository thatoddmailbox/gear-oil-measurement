import nanovnacmd as nv
import lmfit as lm
import numpy as np

from NanoVNASaver.Touchstone import Options, Touchstone

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

def residual(params, x, data=None):
	"""params is the lmfit object. x is the frequency in Hz. For reference: https://lmfit.github.io/lmfit-py/fitting.html
	#based on Ray Kwok & Ji-Fuh Liang. Characterization of High-Q Resonators for Microwave-Filter Applications. IEEE Trans MTT vol.47, p111-114, 1999

	"""
	a = params['a'].value
	b = params['b'].value
	qint = params['qint'].value
	beta = params['beta'].value
	f0 = params['f0'].value
	tau = params['tau'].value
	Omega = x/f0 - f0/x
	S11 = (a + b * 1j) *  ((1-beta)+1j*qint*Omega)/((1+beta)+1j*qint*Omega) * np.exp(1j * (x - f0) * tau)
	# S21 = (a + b * 1j) * (1 - np.exp(1j * phi) * (q / qext) / (1 + 2j * q * ((x - f0) / f0))) * np.exp(1j * (x - f0) * tau)

	if data is None:
		return S11.view(np.float)
	else:
		return (S11.view(np.float)-data.view(np.float))

def ResonatorRegression(frequency=None, smith=None, electricalDelay=None):
	spectrum = smith
	x = frequency
	y = spectrum.view(np.float)
	print(y[0], y[1])

	# create a set of Parameters
	params = lm.Parameters()

	params.add('a', value=np.mean(np.real(spectrum)))
	params.add('b', value=np.mean(np.imag(spectrum)))
	params.add('qint', value=10000, min=0.1)
	#guess at the resonance frequency
	spectrumAbs = np.abs(spectrum)
	frequency_guess = x[np.argmin(spectrumAbs)]
	print("frequency_guess", frequency_guess)
	params.add('f0', value=frequency_guess)
	params.add('tau', value=0)
	params.add('beta', value=1)

	#best Minimizer class documentation i've found so far
	#https://github.com/lmfit/lmfit-py/blob/master/doc/fitting.rst#id81
	mini = lm.Minimizer(residual,params,fcn_args=(x, y))
	return mini

filenames = []

# filenames.append("./data_filings/can3_oil_aluminum_fine.s1p")
# for i in range(1, 13):
# 	filenames.append("./data_filings/can3_oil_aluminum_" + str(i) + "_fine.s1p")
for i in range(0, 8):
	for j in range(0, 2):
		if i == 1 and j == 1:
			continue
		filenames.append("./data_water/water_" + str(i) + "_"+ str(j) + ".s1p")

results = []
for filename in filenames:
	ts = Touchstone(filename)
	ts.load()

	f = np.array([d.freq for d in ts.s11data])
	s11 = np.array([d.z for d in ts.s11data])

	if filename in water_ranges:
		resonance_band_from, resonance_band_to = water_ranges[filename]
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
		f = f[index_from:index_to+1]
		s11 = s11[index_from:index_to+1]

	#running the code
	mini = ResonatorRegression(frequency=f, smith=s11, electricalDelay = None)
	result = mini.minimize()

	f0 = result.params['f0'].value
	a = result.params['a'].value
	b = result.params['b'].value

	qint = result.params['qint'].value
	beta = result.params['beta'].value
	q = qint/(1+beta)
	qext = 1/(1/q-1/qint)

	results.append((f0, a, b, qint, beta, q, qext))

for i, result in enumerate(results):
	print(i, result)