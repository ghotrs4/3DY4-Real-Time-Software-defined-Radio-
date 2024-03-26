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
#include "threadSafeQ.cpp"
#include <iostream>
#include <thread>
#include <cstdio>
#include <queue>
#include <mutex>

using namespace std;

struct RFState {
	std::vector<float> i_state_rf;
	std::vector<float> q_state_rf;
	float prev_I = 0;
	float prev_Q = 0;
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

void frontend(const float rf_decim, const std::vector<float> &rf_coeff, RFState &rf_states, const std::vector<float> &iq_data, threadSafeQ &q)
{
	std::vector<float> i_samples(iq_data.size()/2);
	std::vector<float> q_samples(iq_data.size()/2);

	// separate i and q samples
	for(int i=0;i<iq_data.size();i+=2){
		i_samples[i/2]=iq_data[i];
		q_samples[i/2]=iq_data[i+1];
	}

	std::vector<float> i_downsampled;
	std::vector<float> q_downsampled;
	std::vector<float> fm_demod;

	// downsample and filter i and q data
	downsampleBlockConvolveFIR(rf_decim, i_downsampled, i_samples, rf_coeff, rf_states.i_state_rf);
	downsampleBlockConvolveFIR(rf_decim, q_downsampled, q_samples, rf_coeff, rf_states.q_state_rf);
	
	// use i and q data to obtain demodulated data
	fmDemodArctan(i_downsampled, q_downsampled, rf_states.prev_I, rf_states.prev_Q, fm_demod);

	q.enqueue(fm_demod);
}

void backend(const float audio_Fs, const int audio_decim, const int audio_upsample, const AudioFilter &audio_filters, AudioState &audio_states, PLLState &pll_states, std::vector<float> &audio_block, std::vector<float> &stereo_left, std::vector<float> &stereo_right, threadSafeQ &q)
{
	std::vector<float> pilot_filtered;
	std::vector<float> stereo_filtered;

	// PLL variables
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

	//----------start mono------------

	std::vector<float> dequeuedFMDemod;
    dequeuedFMDemod = q.dequeue();

	delayBlock(dequeuedFMDemod, audio_states.mono_delay_state, mono_delay);

	resampleBlockConvolveFIR(audio_upsample, audio_decim, audio_block, mono_delay, audio_filters.audio_coeff, audio_states.state_audio);
	
	//----------start stereo------------

	blockConvolveFIR(pilot_filtered, dequeuedFMDemod, audio_filters.pilot_coeff, audio_states.pilot_state);
	blockConvolveFIR(stereo_filtered, dequeuedFMDemod, audio_filters.stereo_coeff, audio_states.stereo_state);
	
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


int main(int argc, char* argv[])
{
	int mode = 0;
	bool mono = true;
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

	if (argc == 3) {
		mode = atoi(argv[1]);
		std::string channel = argv[2];
		mono = channel == "stereo" ? false: true;
		if (mode > 3) {
			std::cerr << "Wrong mode: " << mode << std::endl;
			exit(1);
		} else if (channel != "mono" && channel != "stereo") {
			std::cerr << "Wrong parameter: " << channel << ", must be mono or stereo" << std::endl;
		}
	} else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode> <mono/stereo>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}

	std::cerr << "Operating in mode " << mode << (mono ? " mono" : " stereo") << std::endl;

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

			block_size = 1024* audio_decim * rf_decim *2;
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
			rf_decim = 5;

			audio_Fs = 384e3;
			audio_decim = 1280;
			audio_upsample = 147;
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
	std::vector<float> audio_data;
	std::vector<float> stereo_data_left;
	std::vector<float> stereo_data_right;
	/* ------------------------- */
	
	std::vector<float> audio_data_block;
	std::vector<float> stereo_data_left_block;
	std::vector<float> stereo_data_right_block;

	std::vector<float> iq_data(block_size);
	std::vector<float> processed_data(block_size/(audio_decim*rf_decim*2)*audio_upsample * (mono ? 1 : 2));
	std::vector<short int> final_data(block_size/(audio_decim*rf_decim*2)*audio_upsample * (mono ? 1 : 2));

	for (unsigned int block_id=0; ; block_id++) {
		std::cerr << "Block number " << block_id << std::endl;

		readStdinBlockData(block_size, block_id, iq_data);
		if ((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			exit(1);
		}
		
		// do front-end stuff (get fm_demod)
		std::thread frontendThread = std::thread(frontend, std::ref(rf_decim), std::ref(rf_coeff), std::ref(rf_states), std::ref(iq_data), std::ref(q));
        
		// do back-end stuff
        std::thread backendThread = std::thread(backend, std::ref(audio_Fs), std::ref(audio_decim), std::ref(audio_upsample), std::ref(audio_filters), std::ref(audio_states), std::ref(pll_states), std::ref(audio_data_block), std::ref(stereo_data_left_block), std::ref(stereo_data_right_block), std::ref(q));
        
        frontendThread.join();
        backendThread.join();

		if (mono) {
			processed_data.assign(audio_data_block.begin(), audio_data_block.end());
		} else {
			interleave(stereo_data_left_block, stereo_data_right_block, processed_data);
		}
		
		for (unsigned int k=0; k<processed_data.size(); k++) {
			if (std::isnan(processed_data[k])) final_data[k] = 0;
			else final_data[k] = static_cast<short int>(processed_data[k] * 16384);
		}
		fwrite(&final_data[0], sizeof(short int), final_data.size(), stdout);
	}
	
	// std::cout <<"size of output: "<<audio_data.size()<<std::endl;
	// std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	// const std::string out_fname = "../data/float32samples.bin";
	// const std::string out_fname_stereo = "../data/float32samplesStereo.bin";

	// writeBinData(out_fname,audio_data);
	// write_audio_data(out_fname_stereo, stereo_data_left, stereo_data_right);

	return 0;
}
