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
from fmSupportLib import fmDemodArctan, fmPlotPSD, plotSamples
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
	ncoOut[0] = state.ncoState
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
		ncoOut[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)
	state.ncoState = ncoOut[len(pllIn)]
	ncoOut=ncoOut[:-1]
	return ncoOut
def delayBlock(input_block, state_block):
	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]

	return output_block, state_block
def pointwiseMultiply(input1, input2):
	output = np.zeros(len(input1))
	for i in range(len(input1)):
		output[i] = input1[i]*input2[i]*2
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
def manchesterEncoder():
	pass
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
	subfig_height = np.array([2, 2, 2]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 25 * RDS_decim * rf_decim * audio_decim * 2
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
		if il_vs_th == 0:
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

		ncoOut = fmPll(pilot_filtered, pll_freq, audio_Fs, ncoScale, phaseAdjust, normBandwidth, pllState)
		stereo_mixed = pointwiseMultiply(ncoOut, stereo_filtered)

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
		ncoOutRDS = fmPll(RDS_Carrier_filtered, pll_freq_RDS, audio_Fs, ncoScale_RDS, phaseAdjust_RDS, normBandwidth_RDS, pllStateRDS)
		RDS_mixed = pointwiseMultiply(ncoOutRDS, RDS_delay)

		# #Outputs the downsampled and low-pass filtered RDS
		RDS_lowpass, RDS_lowpass_state = resampler(RDS_upsample, RDS_decim, RDS_mixed, RDS_lpf_coeff, RDS_lowpass_state)

		# #Generates root raised cosine impulse response
		RRC_Impulse = impulseResponseRootRaisedCosine(RDS_Fs, RDS_taps)
		RRC_Final, RRC_state = convolve(RDS_lowpass, RRC_Impulse, RRC_state)


		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 4 and block_count < 5:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			#fmPlotPSD(ax0, audio_block, (final_Fs)/1e3, subfig_height[0], \
			#		'PSD audio_block (block ' + str(block_count) + ')')
			plotSamples(ax0, RRC_Final, 2, 10, "RDS after Root raised cosine")
			# output binary file name (where samples are written from Python)
			#fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			#fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			#audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
			ax1.clear()
			fmPlotPSD(ax1,stereo_lowpass, (final_Fs)/1e3, subfig_height[1], \
					'PSD stereo_lowpass (block ' + str(block_count) + ')')

			# plot PSD of selected block after downsampling mono audio
			#audio_block = audio_filt[::audio_decim]
			ax2.clear()
			fmPlotPSD(ax2, stereo_right_data, (final_Fs)/1e3, subfig_height[2], \
					'PSD stereo_right_data (block ' + str(block_count) + ')')

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")
			exit()

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
