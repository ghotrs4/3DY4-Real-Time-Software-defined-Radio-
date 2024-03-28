#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD, plotSamples, manchesterEncoded
from fmRRC import impulseResponseRootRaisedCosine
# for take-home add your functions

# rf_Fs = 1.92e6;
# audio_Fs = 128e3;
# rf_decim = 15;
# rf_Fc = 60e3;
# audio_decim = 4;
# audio_upsample = 1;
# audio_taps *= audio_upsample;
# audio_Fc = 16e3;
# block_size =  1024* audio_decim * rf_decim *2;
# in_fname = "../data/1920.raw";

# rf_Fs = 1.92e6
# rf_Fc = 60e3
# rf_taps = 101
# rf_decim = 15

# audio_Fs = 128e3
# audio_decim = 4
# audio_upsample = 1
# audio_taps = 101
# audio_Fc = 16e3

# final_Fs = 32e3

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10

audio_Fs = 240e3
audio_decim = 5
audio_upsample = 1
audio_taps = 101
audio_taps*=audio_upsample
audio_Fc = 16e3

final_Fs = 48e3

#RDS paramaters
RDS_decim = 120
RDS_upsample = 19
RDS_taps = 101
RDS_taps*=RDS_upsample
RDS_Fc = 3e3
sps = 16
RDS_Fs = sps*2375

in_fname = "../data/samples3.raw"

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1

# call this for first five blocks
# input: symbol stream, last symbol from previous block of symbols, num errors for method 1, num errors for method 2
# ouput: last symbol from this block of symbols, updated num errors for method 1, num errors for method 2
def find_pattern(symbol_stream, symbol_state, errors1, errors2):
	for i in range(1, len(symbol_stream), 2):
		current1 = symbol_stream[i]
		prev1 = symbol_stream[i-1]
		current2 = symbol_stream[i-1]
		prev2 = symbol_stream[i-2] if i != 1 else symbol_state

		if (current1 == prev1):
			errors1 += 1
		if (current2 == prev2):
			errors2 += 1

	symbol_state = symbol_stream[-1]

	return symbol_state, errors1, errors2 # if errors1 > errors2, use start = 0; else use start = 1

# call this after first five blocks
# input: symbol stream, last symbol from previous block of symbols, last bit from previous output bitstream, starting index depending on symbol matching method (0 or 1)
# output: decoded bistream, last symbol from this block of symbols, last bit from this Manchester-decoded (not differentially-decoded) stream
def decode(symbol_stream, symbol_state, bit_state, start):
	bit_stream = []
	output = []

	# Manchester and differential decoding
	for i in range(start, len(symbol_stream), 2):
		current_symbol = symbol_stream[i]
		prev_symbol = symbol_stream[i-1] if i != 0 else symbol_state

		if (current_symbol == 0 and prev_symbol == 1):
			bit_stream.append(1)
			output.append(1 if 1 != bit_state else 0)
		elif (current_symbol == 1 and prev_symbol == 0):
			bit_stream.append(0)
			output.append(1 if 0 != bit_state else 0)
		else:
			# what to do in this case?
			bit_stream.append(0)
			print("error")
		
		bit_state = bit_stream[-1]

	symbol_state = symbol_stream[-1]
	
	return output, symbol_state, bit_state

# input: decoded data, index for next element to be added, boolean indicating whether frames are synced, state from previous window
# output: 16-bit window, index of next element to be added to window, stop flag
# data: input symbol stream, window: output bitstream (length 16), index: index of newest element in window (0 initially), stop: boolean to track when to stop
def get_window(data, index, synced, state):
	index += 16 if synced == True else 1
	# wrap around
	if index >= len(data):
		index -= len(data)

	if index < 15:
		window = state[index:] + data[:index+1]

	else:
		window = data[index-15:index+1]

	state = data[len(data)-15:]

	return window, index, state

# input: 16-bit bitstream
# output: packet of 4 frames (each a 16-bit message stream + 10-bit parity checkword), message with parity but without offset
def frame_sync_transmitter(b):
	packet = EmptyObject()
	p = [0] * 10
	o_a = [0,0,1,1,1,1,1,1,0,0]
	o_b = [0,1,1,0,0,1,1,0,0,0]
	o_c = [0,1,0,1,1,0,1,0,0,0]
	o_cp = [1,1,0,1,0,1,0,0,0,0]
	o_d = [0,1,1,0,1,1,0,1,0,0]

	p[0] = (b[1]+b[2]+b[3]+b[4]+b[5]+b[10]+b[11]+b[12]+b[13]+b[14])%2
	p[1] = (b[2]+b[3]+b[4]+b[5]+b[6]+b[11]+b[12]+b[13]+b[14]+b[15])%2
	p[2] = (b[1]+b[2]+b[6]+b[7]+b[10]+b[11]+b[15])%2
	p[3] = (b[0]+b[1]+b[4]+b[5]+b[7]+b[8]+b[10]+b[13]+b[14])%2
	p[4] = (b[0]+b[1]+b[2]+b[5]+b[6]+b[8]+b[9]+b[11]+b[14]+b[15])%2
	p[5] = (b[0]+b[4]+b[5]+b[6]+b[7]+b[9]+b[11]+b[13]+b[14]+b[15])%2
	p[6] = (b[2]+b[3]+b[4]+b[6]+b[7]+b[8]+b[11]+b[13]+b[15])%2
	p[7] = (b[0]+b[1]+b[2]+b[7]+b[8]+b[9]+b[10]+b[11]+b[13])%2
	p[8] = (b[0]+b[1]+b[2]+b[3]+b[8]+b[9]+b[10]+b[11]+b[12]+b[14])%2
	p[9] = (b[0]+b[1]+b[2]+b[13]+b[4]+b[9]+b[10]+b[11]+b[12]+b[13]+b[15])%2

	packet.a = b + xor(p, o_a)
	packet.b = b + xor(p, o_b)
	packet.c = b + xor(p, o_c)
	packet.cp = b + xor(p, o_cp)
	packet.d = b + xor(p, o_d)

	return packet

# input: packet of four 26-bit frames
def frame_sync_receiver(packet, synced, offsetState, numErrors, position):
	print("new packet received")
	msgs = EmptyObject()
	msgs.a = []
	msgs.b = []
	msgs.c = []
	msgs.d = []

	s = [0] * 10 # syndrome
	m = [] # 26 bit message + offset

	for i in range(1, 6):
		if i == 1:
			m = packet.a
		elif i == 2:
			m = packet.b
		elif i == 3:
			m = packet.c
		elif i == 4:
			m = packet.cp
		elif i == 5:
			m = packet.d
		
		#Defines variables after matrix multiplication
		s[0] = (m[0] + m[10] + m[13] + m[14] + m[15] + m[16] + m[17] + m[19] + m[20] + m[23] + m[24] + m[25])%2
		s[1] = (m[1] + m[11] + m[14] + m[15] + m[16] + m[17] + m[18] + m[20] + m[21] + m[24] + m[25])%2
		s[2] = (m[2] + m[10] + m[12] + m[13] + m[14] + m[18] + m[20] + m[21] + m[22] + m[23] + m[24])%2
		s[3] = (m[3] + m[10] + m[11] + m[16] + m[17] + m[20] + m[21] + m[22])%2
		s[4] = (m[4] + m[11] + m[12] + m[17] + m[18] + m[21] + m[22] + m[23])%2
		s[5] = (m[5] + m[10] + m[12] + m[14] + m[15] + m[16] + m[17] + m[18] + m[20] + m[22] + m[25])%2
		s[6] = (m[6] + m[10] + m[11] + m[14] + m[18] + m[20] + m[21] + m[24] + m[25])%2
		s[7] = (m[7] + m[10] + m[11] + m[12] + m[13] + m[14] + m[16] + m[17] + m[20] + m[21] + m[22] + m[23]+ m[24])%2
		s[8] = (m[8] + m[11] + m[12] + m[13] + m[14] + m[15] + m[17] + m[18] + m[21] + m[22] + m[23] + m[24]+ m[25])%2
		s[9] = (m[9] + m[12] + m[13] + m[14] + m[15] + m[16] + m[18] + m[19] + m[22] + m[23] + m[24] + m[25])%2

		#Checking for offset A		
		if i == 1:
			print("generated type A syndrome: ", s)
			if s == [1,1,1,1,0,1,1,0,0,0]:
				synced = True
				msgs.a = m[0:16]
				print("---------------Block type A found! Bit position ", position)
				if offsetState != 'D' and offsetState != '':
					# false positive
					print("---------------Likely false positive with A syndrome starting at position ", position)
					numErrors += 1
			else:
				numErrors += 1
			offsetState = 'A' if synced else ''
		
		#Checking for offset B
		elif i == 2:
			print("generated type B syndrome: ", s)
			if s == [1,1,1,1,0,1,0,1,0,0]:
				synced = True
				msgs.b = m[0:16]
				print("---------------Block type B found! Bit position ", position)
				if offsetState != 'A' and offsetState != '':
					# false positive
					print("---------------Likely false positive with B syndrome starting at position ", position)
					numErrors += 1
			else:
				numErrors += 1
			offsetState = 'B' if synced else ''
		
		#Checking for offset C
		elif i == 3:
			print("generated type C syndrome: ", s)
			if s == [1,0,0,1,0,1,1,1,0,0]:
				synced = True
				msgs.c = m[0:16]
				print("---------------Block type C found! Bit position ", position)
				if offsetState != 'B' and offsetState != '':
					# false positive
					print("---------------Likely false positive with C syndrome starting at position ", position)
					numErrors += 1
				i += 1 # skip Cp check if C found
			else:
				numErrors += 1
			offsetState = 'C' if synced else ''
		
		#Checking for offset C'
		elif i == 4:
			print("generated type Cp syndrome: ", s)
			if s == [1,1,1,1,0,0,1,1,0,0]:
				synced = True
				msgs.c = m[0:16]
				print("---------------Block type Cp found! Bit position ", position)
				if offsetState != 'B' and offsetState != '':
					# false positive
					print("---------------Likely false positive with Cp syndrome starting at position ", position)
					numErrors += 1
			else:
				numErrors += 1
			offsetState = 'Cp' if synced else ''
			position -= 26

		#Checking for offset D
		elif i == 5:
			print("generated type D syndrome: ", s)
			if s == [1,0,0,1,0,1,1,0,0,0]:
				synced = True
				msgs.d = m[0:16]
				print("---------------Block type D found! Bit position ", position)
				if offsetState != 'C' and offsetState != 'Cp' and offsetState != '':
					# false positive
					print("---------------Likely false positive with D syndrome starting at position ", position)
					numErrors += 1
			else:
				numErrors += 1
			offsetState = 'D' if synced else ''
		
		position += 26
	
	if numErrors > 10 or not synced:
		synced = False
		numErrors = 0
		offsetState = ''
		print("!!!SYNCING!!!")
		
	return synced, msgs, offsetState, numErrors, position

# assumes a and b are arrays of equal length
def xor(a, b):
	result = [0] * len(a)
	for i in range(len(a)):
		result[i] = 1 if (a[i] and not b[i]) or (not a[i] and b[i]) else 0
	return result

def fmDemodArctanCustom(I, Q, prev_I=0, prev_Q=0):
	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):
		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		if(I[k]**2+Q[k]**2 ==0):
			fm_demod[k]=0
			continue
		if(k>0):
			fm_demod[k] = ((I[k]*(Q[k]-Q[k-1]))-(Q[k]*(I[k]-I[k-1])))/(I[k]**2+Q[k]**2)
		else:
			fm_demod[k] = ((I[k]*(Q[k]-prev_Q))-(Q[k]*(I[k]-prev_I)))/(I[k]**2+Q[k]**2)

		# save the state of the current phase
		# to compute the next derivative
	prev_I=I[-1]
	prev_Q=Q[-1]

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_I, prev_Q

def convolve(xb, h, i_state):
	yb = np.zeros(len(xb))
	for n in range(len(yb)):
		for k in range(len(h)):
			if(n-k>=0):
				yb[n]+=h[k]*xb[n-k]
			else:
				yb[n]+=h[k]*i_state[len(i_state) + (n-k)]
	state=xb[len(xb)-(len(h)-1):len(xb)]
	return yb,state

def resampler(upFactor, downFactor, xb, h, i_state):
	yb = np.zeros(int(len(xb)/downFactor*upFactor))
	for n in range(0, len(xb)*upFactor, downFactor):
		phase = n%upFactor
		for k in range(phase, len(h), upFactor):
			if(n-k>=0):
				yb[int(n/downFactor)]+=h[k]*xb[int((n-k)/upFactor)]
			else:
				yb[int(n/downFactor)]+=h[k]*i_state[len(i_state) + int((n-k)/upFactor)]
	state=xb[len(xb)-(int(len(h)/upFactor)-1):len(xb)]
	return yb,state
class EmptyObject:
    pass
def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01, state = EmptyObject()):
	Cp = 2.666
	Ci = 3.555

	Kp = (normBandwidth)*Cp
	Ki = (normBandwidth*normBandwidth)*Ci

	ncoOut = np.empty(len(pllIn)+1)
	q_ncoOut = np.empty(len(pllIn)+1)

	ncoOut[0] = state.ncoState
	q_ncoOut[0] = state.q_ncoState
	for k in range(len(pllIn)):
		errorI = pllIn[k] * (+state.feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-state.feedbackQ)  # feedback complex exponential
		if(errorI==0): #touched
			errorD=0
		else:
			errorD = math.atan2(errorQ, errorI)
	
		state.integrator = state.integrator + Ki*errorD
		state.phaseEst = state.phaseEst + Kp*errorD + state.integrator

		# internal oscillator
		state.trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(state.trigOffset) + state.phaseEst
		state.feedbackI = math.cos(trigArg)
		state.feedbackQ = math.sin(trigArg)

		#return in-phase for stereo, both I and Q for RDS...
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		q_ncoOut[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)
	state.ncoState = ncoOut[len(pllIn)]
	state.q_ncoState = q_ncoOut[len(pllIn)]

	return ncoOut[:-1], q_ncoOut[:-1]
def delayBlock(input_block, state_block):
	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]

	return output_block, state_block
def pointwiseMultiply(input1, input2, gain):
	output = np.zeros(len(input1))
	for i in range(len(input1)):
		output[i] = input1[i]*input2[i]*gain
	return output
def pointwiseAdd(input1, input2):
	output = np.zeros(len(input1))
	for i in range(len(input1)):
		output[i] = input1[i]+input2[i]
	return output
def pointwiseSubtract(input1, input2):
	output = np.zeros(len(input1))
	for i in range(len(input1)):
		output[i] = input1[i]-input2[i]
	return output
def squaringNonlinearity(x):
    y = []  # Initialize an empty list to store the squared values
    for i in range(len(x)):
        y.append(x[i] * x[i])  # Square each element and append to y
    return y

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	raw_data = np.fromfile(in_fname, dtype='uint8')
	# for i in range(5):
	# 	print(raw_data[i])
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
	# init PLL variables
	pll_freq=19e3
	ncoScale=2.0
	phaseAdjust=0
	normBandwidth=0.01

	#init PLL states
	pllState = EmptyObject()
	pllState.integrator = 0.0
	pllState.phaseEst = 0.0
	pllState.feedbackI = 1.0
	pllState.feedbackQ = 0.0
	pllState.ncoState = 1.0
	pllState.trigOffset = 0
	ncoOut = np.zeros(1)

	pllState.q_ncoState = 1.0

	#RDS variabls
	pll_freq_RDS = 114e3
	ncoScale_RDS = 0.5
	phaseAdjust_RDS = 0
	normBandwidth_RDS = 0.003

	#init RDS PLL states
	pllStateRDS = EmptyObject()
	pllStateRDS.integrator = 0.0
	pllStateRDS.phaseEst = 0.0
	pllStateRDS.feedbackI = 1.0
	pllStateRDS.feedbackQ = 0.0
	pllStateRDS.ncoState = 1.0
	pllStateRDS.trigOffset = 0
	ncoOutRDS = np.zeros(1)

	pllStateRDS.q_ncoState = 1.0
	q_ncoOutRDS = np.zeros(1)

	#band pass filter coefficients
	lowcut = 22e3/(audio_Fs/2)
	highcut = 54e3/(audio_Fs/2)
	stereo_coeff = signal.firwin(rf_taps, [lowcut, highcut], window=('hann'), pass_zero=False)
	stereo_state = np.zeros(rf_taps-1)
	stereo_filtered = np.zeros(1)

	lowcut = 18.5e3/(audio_Fs/2)
	highcut = 19.5e3/(audio_Fs/2)
	pilot_coeff = signal.firwin(rf_taps, [lowcut, highcut], window=('hann'), pass_zero=False)
	pilot_state = np.zeros(rf_taps-1)
	pilot_filtered = np.zeros(1)

	mono_delay_state = np.zeros(int(rf_taps/2))#touched, should be -1?
	mono_delay = np.zeros(5) #arbitrary, gets resized

	stereo_mixed = np.zeros(1)
	stereo_lowpass = np.zeros(1)
	stereo_lowpass_state = np.zeros(rf_taps-1)

	stereo_left = np.zeros(1)
	stereo_right = np.zeros(1)

	stereo_left_data = np.empty(1)
	stereo_right_data = np.empty(1)

	# #RDS Channel Extraction Variables
	lowcut = 54e3/(audio_Fs/2)
	highcut = 60e3/(audio_Fs/2)
	RDS_coeff = signal.firwin(rf_taps, [lowcut, highcut], window=('hann'), pass_zero=False)
	RDS_state = np.zeros(rf_taps-1)
	RDS_filtered = np.zeros(1)
	

	# #RDS Carrier Extraction Variables
	lowcut = 113.5e3/(audio_Fs/2)
	highcut = 114.5e3/(audio_Fs/2)
	RDS_Carrier_coeff = signal.firwin(rf_taps, [lowcut, highcut], window=('hann'), pass_zero=False)
	RDS_Squared = np.zeros(1)
	RDS_Carrier_state = np.zeros(rf_taps-1)
	RDS_Carrier_filtered = np.zeros(1)

	# #RDS Delay Path Variables
	RDS_delay_state = np.zeros(int(rf_taps/2))
	RDS_delay = np.zeros(5) #arbitrary, gets resized

	# #RDS Modulation Variables
	RDS_mixed = np.zeros(1)
	RDS_lowpass = np.zeros(1)
	RDS_lowpass_state = np.zeros(RDS_taps-1)


	RDS_lpf_coeff = signal.firwin(RDS_taps*RDS_upsample, RDS_Fc/(audio_Fs*RDS_upsample/2))
	RDS_lpf_coeff=RDS_lpf_coeff*RDS_upsample

	# #RRC Variables
	RRC_Final = np.zeros(1)
	RRC_state = np.zeros(RDS_taps-1)

	#quadrature debug
	q_RDS_mixed = np.zeros(1)
	q_RDS_lowpass = np.zeros(1)
	q_RDS_lowpass_state = np.zeros(RDS_taps-1)

	# #RRC Variables
	q_RRC_Final = np.zeros(1)
	q_RRC_state = np.zeros(RDS_taps-1)
	#RDS
	RDS_encoded = np.empty(1)
	QRDS_encoded = np.empty(1)
	m_index = 0.0
	m_found = False
	RDS_symbols = np.array([])

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs*audio_upsample/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		h = np.zeros(audio_taps)
		Norm_cutoff=audio_Fc/(audio_Fs*audio_upsample/2)
		for i in range(audio_taps):
			if (i == (audio_taps-1)/2):
				h[i] = Norm_cutoff
			else:
				param = np.pi * Norm_cutoff * (i - (audio_taps - 1)/2)
				# print("param: " + str(param))
				h[i] = Norm_cutoff * np.sin(param)/param
			h[i] = h[i]*((np.sin(i*np.pi/audio_taps))**2)*audio_upsample #gain
		audio_coeff = h

	# set up the subfigures for plotting
	subfig_height = np.array([2, 2, 4, 2]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size =  sps * RDS_decim * rf_decim * audio_decim * 2 * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0
	# add state as needed for the mono channel filter

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# custom resets
	prev_I = 0
	prev_Q = 0
	state = np.zeros(audio_taps-1)

	window_index = 30
	synced = False
	window_state = []

	symbol_state = 0
	errors1 = 0
	errors2 = 0

	decode_start = 0

	bit_state = 0

	offsetState = ''
	numErrors = 0
	bit_pos = 0

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		print('Processing block ' + str(block_count))
		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 1:
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, prev_I, prev_Q = fmDemodArctanCustom(i_ds, q_ds, prev_I, prev_Q)


		# extract the mono audio data through filtering
		if il_vs_th == 0:
			# to be updated by you during the in-lab session based on lfilter
			# same principle as for i_filt or q_filt (but different arguments)
			# print("audio coeff size: " + str(len(audio_coeff)))
			# print("fm_demod size: " + str(len(fm_demod)))
			# print("state size: " + str(len(state)))
			audio_filt, state = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=state)
		else:
			# to be updated by you for the takehome exercise
			# with your own code for BLOCK convolution
			#audio_filt, state = convolve(fm_demod, audio_coeff, state)
			mono_delay, mono_delay_state = delayBlock(fm_demod, mono_delay_state)
			audio_block, state = resampler(audio_upsample, audio_decim, mono_delay, audio_coeff, state)
		
		pilot_filtered, pilot_state = convolve(fm_demod, pilot_coeff, pilot_state)
		stereo_filtered, stereo_state = convolve(fm_demod, stereo_coeff, stereo_state)

		ncoOut, _ = fmPll(pilot_filtered, pll_freq, audio_Fs, ncoScale, phaseAdjust, normBandwidth, pllState)
		stereo_mixed = pointwiseMultiply(ncoOut, stereo_filtered, 2)

		stereo_lowpass, stereo_lowpass_state = resampler(audio_upsample, audio_decim, stereo_mixed, audio_coeff, stereo_lowpass_state)

		stereo_left = pointwiseAdd(audio_block, stereo_lowpass)
		stereo_right = pointwiseSubtract(audio_block, stereo_lowpass)

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		#audio_data = np.concatenate((audio_data, audio_block))
		stereo_left_data = np.concatenate((stereo_left_data, stereo_left))
		stereo_right_data = np.concatenate((stereo_right_data, stereo_right))
		#

		#-------------RDS---------------
		#RDS Channel Extraction
		RDS_filtered, RDS_state = convolve(fm_demod, RDS_coeff, RDS_state)

		# #RDS Carrier Extraction
		RDS_Squared = squaringNonlinearity(RDS_filtered)
		RDS_Carrier_filtered, RDS_Carrier_state = convolve(RDS_Squared, RDS_Carrier_coeff, RDS_Carrier_state)

		# #RDS Delay Path
		RDS_delay, RDS_delay_state = delayBlock(RDS_filtered, RDS_delay_state)

		# #RDS pll
		ncoOutRDS, q_ncoOutRDS = fmPll(RDS_Carrier_filtered, pll_freq_RDS, audio_Fs, ncoScale_RDS, phaseAdjust_RDS, normBandwidth_RDS, pllStateRDS)
		RDS_mixed = pointwiseMultiply(ncoOutRDS, RDS_delay, 1)

		# #Outputs the downsampled and low-pass filtered RDS
		RDS_lowpass, RDS_lowpass_state = resampler(RDS_upsample, RDS_decim, RDS_mixed, RDS_lpf_coeff, RDS_lowpass_state)

		# #Generates root raised cosine impulse response
		RRC_Impulse = impulseResponseRootRaisedCosine(RDS_Fs, RDS_taps)
		RRC_Final, RRC_state = convolve(RDS_lowpass, RRC_Impulse, RRC_state)

		#quadrature debug path
		q_RDS_mixed = pointwiseMultiply(q_ncoOutRDS, RDS_delay, 1)
		q_RDS_lowpass, q_RDS_lowpass_state = resampler(RDS_upsample, RDS_decim, q_RDS_mixed, RDS_lpf_coeff, q_RDS_lowpass_state)
		q_RRC_Final, q_RRC_state = convolve(q_RDS_lowpass, RRC_Impulse, q_RRC_state)


		RDS_encoded, QRDS_encoded, RDS_symbols, m_index, m_found = manchesterEncoded(RRC_Final, q_RRC_Final, sps, m_index, m_found)

		if (block_count < 5): # find pattern
				symbol_state, errors1, errors2 = find_pattern(RDS_symbols, symbol_state, errors1, errors2)
		else: # decode
			decode_start = 0 if errors1 > errors2 else 1
			decoded_stream, symbol_state, bit_state = decode(RDS_symbols, symbol_state, bit_state, decode_start)

			widx = 0
			while ((synced and widx < int(len(RDS_symbols)/2)-16) or (not synced and widx < int(len(RDS_symbols)/2))):
				window_data, window_index, window_state = get_window(decoded_stream, window_index, synced, window_state)
				widx = window_index

				packet = frame_sync_transmitter(window_data)
				synced, msgs, offsetState, numErrors, bit_pos = frame_sync_receiver(packet, synced, offsetState, numErrors, bit_pos)
				
				print("msgs.a: ", msgs.a)

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		# if block_count >= 40 and block_count < 41:

		# 	# plot PSD of selected block after FM demodulation
		# 	ax0.clear()
		# 	#fmPlotPSD(ax0, audio_block, (final_Fs)/1e3, subfig_height[0], \
		# 	#		'PSD audio_block (block ' + str(block_count) + ')')
		# 	plotSamples(ax0, RRC_Final, 2, 1, "RDS RRC In-phase")
		# 	# output binary file name (where samples are written from Python)
		# 	#fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
		# 	# create binary file where each sample is a 32-bit float
		# 	#fm_demod.astype('float32').tofile(fm_demod_fname)

		# 	# plot PSD of selected block after extracting mono audio
		# 	#audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
		# 	ax1.clear()
		# 	plotSamples(ax1, q_RRC_Final, 2, 1, "RDS RRC Quadrature")
		# 	#fmPlotPSD(ax1,stereo_lowpass, (final_Fs)/1e3, subfig_height[1], \
		# 	#		'PSD stereo_lowpass (block ' + str(block_count) + ')')

		# 	# plot PSD of selected block after downsampling mono audio
		# 	#audio_block = audio_filt[::audio_decim]
		# 	ax2.clear()
		# 	print("Manchester encoded symbols: "+str(len(RDS_encoded)))
		# 	ax2.scatter(RDS_encoded, QRDS_encoded, s=10)
		# 	#fmPlotPSD(ax2, stereo_right_data, (final_Fs)/1e3, subfig_height[2], \
		# 			#'PSD stereo_right_data (block ' + str(block_count) + ')')
		# 	ax3.clear()
		# 	#fmPlotPSD(ax0, audio_block, (final_Fs)/1e3, subfig_height[0], \
		# 	#		'PSD audio_block (block ' + str(block_count) + ')')
		# 	x_vec = np.zeros(len(RDS_symbols))
		# 	for i in range(len(RDS_symbols)):
		# 		x_vec[i]=i
				
		# 	ax3.scatter(x_vec, RDS_symbols, s=10)
		# 	# save figure to file
		# 	fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")
			# exit()

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	#wavfile.write(out_fname, int(final_Fs), np.int16((audio_data/2)*32767))
	stereo_data = np.vstack((stereo_left_data, stereo_right_data)).T
	wavfile.write(out_fname, int(final_Fs), np.int16((stereo_data/2)*32767))#touched
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	plt.show()
