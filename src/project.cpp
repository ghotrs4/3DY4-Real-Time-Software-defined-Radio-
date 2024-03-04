/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

using namespace std;

void mono(const int mode)
{
	//for debug, add mode functionality
	float rf_Fs = 2.4e6;
	float rf_Fc = 100e3;
	unsigned short int rf_taps = 101;
	float rf_decim = mode == 0 ? 10 : 5;//add modes

	float audio_Fs = 48e3;
	float audio_decim = mode == 0 ? 5 : 8;//add modes
	unsigned short int audio_taps = 101;
	float audio_Fc = 16e3;

	std::vector<float> audio_data; //output audio sample vectpr

	const std::string in_fname = "../data/iq_samples.raw";
	std::vector<uint8_t> raw_data;
	readRawData(in_fname, raw_data);
	std::vector<double> iq_data;
	convertRaw(raw_data, iq_data);
	cout << iq_data[0] << endl;

	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	cout<<"rf coefficients generated!"<<endl;//debug

	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff);
	cout<<"impulseResponseLPF generated!"<<endl;//debug
	
	std::vector<float> y;
	std::vector<float> state;
	int position = 0;
	int block_size = 1000;
	std::vector<float> i_samples;
	std::vector<float> q_samples;
	float prev_I=0;
	float prev_Q=0;
	std::vector<float> fm_demod;

	for(int i=0;i<iq_data.size();i+=2){
		i_samples.push_back(iq_data[i]);
		q_samples.push_back(iq_data[i+1]);
	}
	cout<<"here 1"<<endl;//debug

	while ((position+1) * block_size < iq_data.size()) {
		std::vector<float> i_block;
		std::vector<float> q_block;

		blockConvolveFIR(i_block, i_samples, rf_coeff, state, position, block_size);
		blockConvolveFIR(q_block, q_samples, rf_coeff, state, position, block_size);
		cout<<"convolve 1"<<endl;//debug

		// downsample
		
		std::vector<float> i_downsampled;
		std::vector<float> q_downsampled;
		downsample(i_block, rf_decim, i_downsampled);
		downsample(q_block, rf_decim, q_downsampled);
		cout<<"downsample 1"<<endl;//debug

		fmDemodArctan(i_downsampled, q_downsampled, prev_I, prev_Q, fm_demod);
		cout<<"demod"<<endl;//debug

		std::vector<float> audio_filt;
		cout<<"fm_demod size:"<<fm_demod.size()<<std::endl;
		cout<<"audio_coeff size:"<<audio_coeff.size()<<std::endl;
		blockConvolveFIR(audio_filt, fm_demod, audio_coeff, state, position, block_size);
		cout<<"convolve 2"<<endl;//debug

		position += block_size;

		std::vector<float> audio_block;
		downsample(audio_filt, audio_decim, audio_block);
		cout<<"downsample 2"<<endl;//debug

		audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());
		cout<<"insert"<<endl;//debug
	}
}

int main()
{
	int mode = 0;
	mono(mode);
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
	std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp
	DFT(slice_data, Xf);

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp
	computeVectorMagnitude(Xf, Xmag);

	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//
	std::vector<float> freq;
	std::vector<float> psd_est;
	estimatePSD(bin_data, NFFT, 240, freq, psd_est);
	logVector("demod_psd", freq, psd_est); // log only positive freq

	// if you wish to write some binary files, see below example
	//
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);
	//
	// output files can be imported, for example, in Python
	// for additional analysis or alternative forms of visualization

	// naturally, you can comment the line below once you are comfortable to run GNU plot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
