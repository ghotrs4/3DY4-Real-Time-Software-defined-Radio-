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

struct RFState {
	std::vector<float> i_state_rf;
	std::vector<float> q_state_rf;
};

struct AudioState {
	std::vector<float> state_audio;
	std::vector<float> stereo_lowpass_state;
	std::vector<float> pilot_state;
	std::vector<float> stereo_state;
	std::vector<float> mono_delay_state;
};

struct AudioFilter {
	std::vector<float> audio_coeff;
	std::vector<float> pilot_coeff;
	std::vector<float> stereo_coeff;
};

struct PLLState {
    float feedbackI = 1.0;
    float feedbackQ = 0.0;
    float integrator = 0;
    float phaseEst = 0;
    float trigOffset = 0;
    float nco_state = 1.0;
};

void frontend(const int mode, const float rf_decim, const std::vector<float> &rf_coeff, RFState &rf_states, const std::vector<uint8_t> &raw_data, std::vector<float> &fm_demod)
{
	// obtain iq data
	std::vector<float> iq_data;
	convertRaw(raw_data, iq_data);

	std::vector<float> i_samples;
	std::vector<float> q_samples;
	float prev_I=0;
	float prev_Q=0;

	// separate i and q samples
	for(int i=0;i<iq_data.size();i+=2){
		i_samples.push_back(iq_data[i]);
		q_samples.push_back(iq_data[i+1]);
	}

	std::vector<float> i_downsampled;
	std::vector<float> q_downsampled;

	// downsample and filter i and q data
	downsampleBlockConvolveFIR(rf_decim, i_downsampled, i_samples, rf_coeff, rf_states.i_state_rf);
	downsampleBlockConvolveFIR(rf_decim, q_downsampled, q_samples, rf_coeff, rf_states.q_state_rf);

	// use i and q data to obtain demodulated data
	fmDemodArctan(i_downsampled, q_downsampled, prev_I, prev_Q, fm_demod);
}

void backend(const int mode, const float audio_Fs, const int audio_decim, const int audio_upsample, const AudioFilter &audio_filters, AudioState &audio_states, PLLState &pll_states, const std::vector<float> &fm_demod, std::vector<float> &audio_block, std::vector<float> &stereo_left, std::vector<float> &stereo_right)
{
	std::vector<float> pilot_filtered;
	std::vector<float> stereo_filtered;

	float pll_freq=19e3;
	float ncoScale=2.0;
	float phaseAdjust=0;
	float normBandwidth=0.01;
	
	std::vector<float> ncoOut;

	std::vector<float> stereo_mixed;
	std::vector<float> stereo_lowpass;

	std::vector<float> mono_delay;
	
	std::vector<float> stereo_block;
	//debug vectors
	// std::vector<float> vector_index;
	// std::vector<float> vector_data;

	delayBlock(fm_demod, audio_states.mono_delay_state, mono_delay);

	resampleBlockConvolveFIR(audio_upsample, audio_decim, audio_block, mono_delay, audio_filters.audio_coeff, audio_states.state_audio);
	
	//----------start stereo------------

	blockConvolveFIR(pilot_filtered, fm_demod, audio_filters.pilot_coeff, audio_states.pilot_state);
	blockConvolveFIR(stereo_filtered, fm_demod, audio_filters.stereo_coeff, audio_states.stereo_state);
	
	fmPLL(pilot_filtered, pll_freq, audio_Fs, ncoScale, phaseAdjust, normBandwidth, ncoOut, pll_states.feedbackI, pll_states.feedbackQ, pll_states.integrator, pll_states.phaseEst, pll_states.trigOffset, pll_states.nco_state);

	//mix carrier + stereo signal
	pointwiseMultiply(ncoOut, stereo_filtered, stereo_mixed);

	//downsample and convolve to achieve desired output sample rate (ie 48k for mode 0)
	resampleBlockConvolveFIR(audio_upsample, audio_decim, stereo_lowpass, stereo_mixed, audio_filters.audio_coeff, audio_states.stereo_lowpass_state);

	//------for plotting don't delete-------
	// if (position == block_size*32) {
	// 	float nfft = 2048;
	// 	std::vector<float> slice_data = \
	// 	std::vector<float>(fm_demod.begin(), fm_demod.begin() + nfft);

	// 	std::vector<std::complex<float>> Xf;
	// 	std::vector<float> Xmag;

	// 	DFT(slice_data, Xf);
	// 	computeVectorMagnitude(Xf, Xmag);
	// 	vector_index.clear();
	// 	genIndexVector(vector_index, Xmag.size());
	// 	for(int i=0;i<vector_index.size();i++){
	// 		vector_index[i]=vector_index[i]*(audio_Fs)/nfft;//change audio_Fs to the correct freqeuncy
	// 	}
	// 	logVector("ncoOut1", vector_index, Xmag); // log only positive freq
	// }

	pointwiseAdd(audio_block, stereo_lowpass, stereo_left);
	pointwiseSubtract(audio_block, stereo_lowpass, stereo_right);
}


int main()
{
	int mode = 3;
	short int num_taps = 101;
	std::string in_fname;
	int block_size;
	// rf params
	float rf_Fs;
	float rf_Fc = 100e3;
	float rf_decim;
	// audio params
	float audio_Fs;
	float audio_Fc = 16e3;
	short int audio_taps;
	float audio_decim;
	float audio_upsample;

	switch(mode) {
		case 0: //output Fs = 48k
			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 240e3;
			audio_decim = 5;
			audio_upsample = 1;
			audio_taps = num_taps;

			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/stereo_l0_r9.raw";
			break;
		case 1://output Fs = 36k
			rf_Fs = 1.44e6;
			rf_decim = 5;

			audio_Fs = 288e3;
			audio_decim = 8;
			audio_upsample = 1;
			audio_taps = num_taps;

			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/1440.raw";
			break;
		case 2://output Fs = 44.1k
			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 240e3;
			audio_decim = 800;
			audio_upsample = 147;
			audio_taps = num_taps * audio_upsample;

			block_size = 10 * audio_decim * rf_decim *2;
			in_fname = "../data/stereo_l0_r9.raw";
			break;
		case 3://output Fs = 44.1k
			rf_Fs = 1.92e6;
			rf_decim = 15;

			audio_Fs = 128e3;
			audio_decim = 1280;
			audio_upsample = 441;
			audio_taps = num_taps * audio_upsample;

			block_size =  10* audio_decim * rf_decim *2	;
			in_fname = "../data/1920.raw";
			break;
		default:
			rf_Fs = 2.4e6;
			rf_decim = 10;

			audio_Fs = 240e3;
			audio_decim = 5;
			audio_upsample = 1;
			audio_taps = num_taps;

			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/2400.raw";
	}

	RFState rf_states;

	rf_states.i_state_rf.resize(num_taps - 1, 0.0);
	rf_states.q_state_rf.resize(num_taps - 1, 0.0);

	AudioState audio_states;

	audio_states.state_audio.resize((num_taps - 1), 0.0);
	audio_states.pilot_state.resize(num_taps-1, 0.0);
	audio_states.stereo_state.resize(num_taps-1, 0.0);
	audio_states.stereo_lowpass_state.resize(num_taps-1, 0.0);
	audio_states.mono_delay_state.resize((num_taps)/2, 0.0);

	AudioFilter audio_filters;

	PLLState pll_states;

	// read raw data (change later for real-time processing)
	std::vector<uint8_t> raw_data;
	readRawData(in_fname, raw_data);

	// obtain filter coefficients for rf low-pass filter
	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, num_taps, rf_coeff, 1);

	// obtain filter coefficients for audio low-pass filter
	impulseResponseLPF(audio_Fs*audio_upsample, audio_Fc, audio_taps, audio_filters.audio_coeff, audio_upsample);

	// obtain filter coefficients for stereo band-pass filters
	float pilot_Fb=18.5e3;
	float pilot_Fe=19.5e3;
	float stereo_Fb=22e3;
	float stereo_Fe=54e3;
	impulseResponseBPF(audio_Fs, pilot_Fb, pilot_Fe, num_taps, audio_filters.pilot_coeff, 1);
	impulseResponseBPF(audio_Fs, stereo_Fb, stereo_Fe, num_taps, audio_filters.stereo_coeff, 1);

	/* ------- temporary (?) ------- */
	int position = 0;
	std::vector<float> fm_demod;
	std::vector<float> audio_data;
	std::vector<float> stereo_data_left;
	std::vector<float> stereo_data_right;
	/* ------------------------- */
	
	std::vector<float> fm_demod_block;
	std::vector<float> audio_data_block;
	std::vector<float> stereo_data_left_block;
	std::vector<float> stereo_data_right_block;
	std::vector<uint8_t> raw_data_block;

	while (position+block_size < raw_data.size()) {
		std::cout<<"block number: "<<position/block_size<<std::endl;

		// get raw data block
		raw_data_block = std::vector<uint8_t>(raw_data.begin() + position, raw_data.begin() + position + block_size);

		// do front-end stuff (get fm_demod)
		frontend(mode, rf_decim, rf_coeff, rf_states, raw_data_block, fm_demod_block);

		// do back-end stuff
		backend(mode, audio_Fs, audio_decim, audio_upsample, audio_filters, audio_states, pll_states, fm_demod_block, audio_data_block, stereo_data_left_block, stereo_data_right_block);

		audio_data.insert(audio_data.end(), audio_data_block.begin(), audio_data_block.end());
		stereo_data_left.insert(stereo_data_left.end(), stereo_data_left_block.begin(), stereo_data_left_block.end());
		stereo_data_right.insert(stereo_data_right.end(), stereo_data_right_block.begin(), stereo_data_right_block.end());

		position += block_size;
	}
	
	std::cout <<"size of output: "<<audio_data.size()<<std::endl;
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	const std::string out_fname = "../data/float32samples.bin";
	const std::string out_fname_stereo = "../data/float32samplesStereo.bin";

	writeBinData(out_fname,audio_data);
	write_audio_data(out_fname_stereo, stereo_data_left, stereo_data_right);

	return 0;
}