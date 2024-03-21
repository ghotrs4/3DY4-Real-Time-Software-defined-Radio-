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

void mono(const int mode,std::vector<float>& audio_data, std::vector<float>& stereo_data_left,std::vector<float>& stereo_data_right)
{
	float rf_Fs;
	float rf_Fc = 100e3;
	unsigned short int rf_taps = 101;
	float rf_decim;

	float audio_Fs;
	int audio_decim;
	int audio_upsample;
	unsigned short int audio_taps = 101;
	float audio_Fc;

	switch(mode) {
		case 0: //output Fs = 48k
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 16e3;
			break;
		case 1://output Fs = 36k
			rf_Fs = 1.44e6;
			audio_Fs = 288e3;
			rf_decim = 5;
			audio_decim = 8;
			audio_upsample = 1;
			audio_Fc = 16e3;
			break;
		case 2://output Fs = 44.1k
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 800;
			audio_upsample = 147;
			audio_taps *= audio_upsample;
			audio_Fs *= audio_upsample;
			audio_Fc = 16e3;
			break;
		case 3://output Fs = 44.1k
			rf_Fs = 1.92e6;
			audio_Fs = 128e3;
			rf_decim = 15;
			audio_decim = 1280;
			audio_upsample = 441;
			audio_taps *= audio_upsample;
			audio_Fs *= audio_upsample;
			audio_Fc = 16e3;
			break;
		default:
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 24e3;
	}

	const std::string in_fname = "../data/iq_samples.raw";
	std::vector<uint8_t> raw_data;
	readRawData(in_fname, raw_data);
	std::vector<float> iq_data;
	convertRaw(raw_data, iq_data);
	cout <<"size of input: "<<raw_data.size()<<endl;
	cout <<"size of iq_data: "<<iq_data.size()<<endl;

	std::vector<float> rf_coeff;
	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff, 1);

	std::vector<float> audio_coeff;
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff, audio_upsample);
	
	std::vector<float> y;

	std::vector<float> i_state_rf;
	std::vector<float> q_state_rf;
	std::vector<float> state_audio;

	i_state_rf.resize(rf_coeff.size() - 1, 0.0);
	q_state_rf.resize(rf_coeff.size() - 1, 0.0);
	state_audio.resize((audio_coeff.size()/audio_upsample - 1), 0.0);

	int position = 0;
	int block_size = 1024 * rf_decim * audio_decim * 2/audio_upsample;
	std::vector<float> i_samples;
	std::vector<float> q_samples;
	float prev_I=0;
	float prev_Q=0;
	std::vector<float> fm_demod;

	for(int i=0;i<iq_data.size();i+=2){
		i_samples.push_back(iq_data[i]);
		q_samples.push_back(iq_data[i+1]);
	}
	cout <<"size of i_samples: "<<i_samples.size()<<endl;
	cout <<"size of q_samples: "<<q_samples.size()<<endl;

	std::vector<float> i_downsampled;
	std::vector<float> q_downsampled;

	std::vector<float> audio_filt;
	std::vector<float> audio_block;

	//stereo variables
	std::vector<float> pilot_coeff;
	std::vector<float> stereo_coeff;

	std::vector<float> pilot_filtered;
	std::vector<float> stereo_filtered;
	float pilot_Fb=18.5e3;
	float pilot_Fe=19.5e3;
	float stereo_Fb=22e3;
	float stereo_Fe=54e3;
	unsigned short int num_taps_stereo=101;

	float pll_freq=19e3;
	float ncoScale=1.0;
	float phaseAdjust=0;
	float normBandwidth=0.01;
	float feedbackI=1.0;
	float feedbackQ=0.0;
	float integrator=0;
	float phaseEst=0;
	float trigOffset=0;
	std::vector<float> ncoOut;

	std::vector<float> stereo_mixed;
	std::vector<float> stereo_lowpass;
	std::vector<float> stereo_lowpass_state;
	
	impulseResponseBPF(audio_Fs,pilot_Fb, pilot_Fe, num_taps_stereo, pilot_coeff, 1);//gain of 1
	impulseResponseBPF(audio_Fs,stereo_Fb, stereo_Fe, num_taps_stereo, stereo_coeff, 1);//gain of 1

	std::vector<float> pilot_state;
	std::vector<float> stereo_state;
	std::vector<float> mono_delay_state;
	
	pilot_state.resize(pilot_coeff.size()-1, 0.0);
	stereo_state.resize(stereo_coeff.size()-1, 0.0);
	stereo_lowpass_state.resize(audio_coeff.size()-1, 0.0);
	mono_delay_state.resize(audio_coeff.size()-1, 0.0);

	std::vector<float> mono_delay;
	
	std::vector<float> stereo_left;
	std::vector<float> stereo_right;
	std::vector<float> stereo_block;

	while (position+block_size<iq_data.size()) {
		cout<<"block number: "<<position/block_size<<endl;

		downsampleBlockConvolveFIR(rf_decim, i_downsampled, i_samples, rf_coeff, i_state_rf, position/2, block_size/2);
		downsampleBlockConvolveFIR(rf_decim, q_downsampled, q_samples, rf_coeff, q_state_rf, position/2, block_size/2);

		fmDemodArctan(i_downsampled, q_downsampled, prev_I, prev_Q, fm_demod);

		//begin mono
		resampleBlockConvolveFIR(audio_upsample, audio_decim, audio_block, fm_demod, audio_coeff, state_audio, 0, fm_demod.size());
		delayBlock(audio_block, mono_delay_state, mono_delay);
		//----------start stereo------------
		
		//isolate pilot and stereo channel
		
		blockConvolveFIR(pilot_filtered, fm_demod, pilot_coeff, pilot_state, 0, fm_demod.size());
		blockConvolveFIR(stereo_filtered, fm_demod, stereo_coeff, stereo_state, 0, fm_demod.size());

		//PLL + NCO to recover carrier (ie pilot tone phase shift to 38kHz)
		fmPLL(pilot_filtered, pll_freq, audio_Fs, ncoScale, phaseAdjust, normBandwidth, ncoOut, feedbackI, feedbackQ, integrator, phaseEst, trigOffset);
		// for (int i = 0; i < 5; i++) {
		// 	cout << "pilot_filtered"
		// }

		//mix carrier + stereo signal
		pointwiseMultiply(ncoOut, stereo_filtered, stereo_mixed);//touched
		// for(int i=0;i<5;i++){
		// 	cout<<"pw multiply output: "<<stereo_mixed[i]<<endl;
		// 	cout<<"nco stereo output: "<<ncoOut[i]<<endl;
		// 	cout<<"stereo_filtered: "<<stereo_filtered[i]<<endl;
		// 	cout<<"pilot_filtered: "<<pilot_filtered[i]<<endl;
		// }
		
		//downsample and convolve to achieve desired output sample rate (ie 48k for mode 0)
		downsampleBlockConvolveFIR(audio_decim, stereo_lowpass, stereo_mixed, audio_coeff, stereo_lowpass_state, 0, stereo_mixed.size());

		//output stereo signal
		pointwiseAdd(stereo_lowpass, mono_delay, stereo_left);
		pointwiseSubtract(stereo_lowpass, mono_delay, stereo_right);
		//prepare new block
		
		//cout<<"position+block_size: "<<position+block_size<<endl;
		//cout<<"iq_data.size(): "<<iq_data.size()<<endl;
		if (position > 0) {//output mono and stereo
			audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());
			stereo_data_left.insert(stereo_data_left.end(), stereo_left.begin(), stereo_left.end());
			stereo_data_right.insert(stereo_data_right.end(), stereo_right.begin(), stereo_right.end());
		}
		position += block_size;
		// cout<<stereo_data_left
	}
}

int main()
{
	int mode = 0;
	std::vector<float> audio_data; //output audio sample vector
	std::vector<float> stereo_data_left;
	std::vector<float> stereo_data_right;

	mono(mode, audio_data, stereo_data_left, stereo_data_right);

	const std::string out_fname = "../data/float32samples.bin";
	const std::string out_fname_stereo = "../data/float32samplesStereo.bin";

	writeBinData(out_fname,audio_data);
	write_audio_data(out_fname_stereo, stereo_data_left, stereo_data_right);
	return 0;
}