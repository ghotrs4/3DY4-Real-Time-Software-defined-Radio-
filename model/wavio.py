# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
from scipy.io import wavfile
import sys

#
# since manipulating .wav files is not the objective of the SDR project and
# we are using them solely for "assessing" the outcome of the DSP tasks
# while troubleshooting, we will avoid processing any .wav files in C++,
# mainly because of the error prone nature of handling .wav file headers
#
# for the reason above, the Python script below can be used to parse/format
# .wav files to/from binary files where the sample representation is known
# (or better said agreed on) by both the Python script and the C++ program
#
# .wav files should be opened only in this Python script and samples written
# in binary (e.g., assuming 32-bit floating point for this example) should be
# read by the C++ program in binary format (raw data, no headers); subsequently,
# the C++ program should output the processed data also in binary formart,
# which can be read back by this Python script to be formatted properly with a
# a header into a .wav file that can then be used on a third part audio player
#

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\tmono:  process mono audio')
	print('\tstereo: process stereo audio')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2 or (sys.argv[1] != 'mono' and sys.argv[1] != 'stereo'):
		cli_error_msg()

	# parse an audio file
	# audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	# print(' Audio sample rate = {0:f} \
	# 	\n Number of channels = {1:d} \
	# 	\n Numbef of samples = {2:d}' \
	# 	.format(audio_Fs, audio_data.ndim, len(audio_data)))

	# output binary file name (where samples are written from Python)
	# out_fname = "../data/float32samples.bin"
	# # dump audio data in a binary file where each sample is a 32-bit float
	# audio_data.astype('float32').tofile(out_fname)
	# print(" Written binary data to \"" + out_fname + "\" in float32 format")

	# input binary file name (from where samples are read into Python)
	# the default is JUST a SELF-CHECK; of course, change filenames as needed
	in_fname = "../data/float32samples.bin" if sys.argv[1] == 'mono' else "../data/float32samplesStereo.bin"
	# in_fname = "../data/float32filtered.bin"
	# read data from a binary file (assuming 32-bit floats)
	float_data = np.fromfile(in_fname, dtype='float32')
	# for i in range(25):
	# 	print(float_data[i])
	print(" Read binary data from \"" + in_fname + "\" in float32 format")

	if (sys.argv[1] == 'stereo'):
		# we assume below there are two audio channels where data is
		# interleaved, i.e., left channel sample, right channel sample, ...
		# for mono .wav files the reshaping below is unnecessary
		reshaped_data = np.reshape(float_data, (-1, 2))

		# # self-check if the read and write are working correctly
		# # not needed while working with data generated from C++
		# print(" Are the two data sets identical ? " +
		# 		str(np.array_equal(float_data,reshaped_data)))

	wavfile.write("../data/audio_processed.wav", \
				int(483), \
				np.int16((float_data)) if not sys.argv[1] == 'stereo' \
				else np.int16((reshaped_data)))

	# note: we can also dump audio data in other formats, if needed
	# audio_data.astype('int16').tofile('int16samples.bin')
