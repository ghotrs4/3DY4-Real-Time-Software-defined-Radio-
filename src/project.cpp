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

void mono(const int mode,std::vector<float>& audio_data)
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
	int block_size;
	std::string in_fname;

	switch(mode) {
		case 0: //output Fs = 48k
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 16e3;
			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/2400.raw";
			break;
		case 1://output Fs = 36k
			rf_Fs = 1.44e6;
			audio_Fs = 288e3;
			rf_decim = 5;
			audio_decim = 8;
			audio_upsample = 1;
			audio_Fc = 18e3;
			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/1440.raw";
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
			block_size = 10 * audio_decim * rf_decim *2;
			in_fname = "../data/2400.raw";
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
			block_size = 10 * audio_decim * rf_decim *2;
			in_fname = "../data/1920.raw";
			break;
		default:
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 24e3;
			block_size = 1024 * audio_decim * rf_decim *2;
			in_fname = "../data/2400.raw";
	}

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
	//2 if for i and q samples, 10 accounts for minimum block sizex
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
	float pilot_Fb=18.5e3;
	float pilot_Fe=19.5e3;
	float stereo_Fb=22e3;
	float stereo_Fe=54e3;
	unsigned short int num_taps_stereo=101;

	impulseResponseBPF(audio_Fs,pilot_Fb, pilot_Fe, num_taps_stereo, pilot_coeff, 1);//gain of 1
	impulseResponseBPF(audio_Fs,stereo_Fb, stereo_Fe, num_taps_stereo, stereo_coeff, 1);//gain of 1


	while (position+block_size<iq_data.size()) {//touched
		cout<<"block number: "<<position/block_size<<endl;

		downsampleBlockConvolveFIR(rf_decim, i_downsampled, i_samples, rf_coeff, i_state_rf, position/2, block_size/2);
		downsampleBlockConvolveFIR(rf_decim, q_downsampled, q_samples, rf_coeff, q_state_rf, position/2, block_size/2);

		fmDemodArctan(i_downsampled, q_downsampled, prev_I, prev_Q, fm_demod);

		resampleBlockConvolveFIR(audio_upsample, audio_decim, audio_block, fm_demod, audio_coeff, state_audio, 0, fm_demod.size());
		if (position != 0) {
			audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());
		}

		//----------start stereo------------

		
		// float feedbackI = 1.0;
		// float feedbackQ = 0.0;
		// std::vector<float> PLLin;
		// std::vector<float> ncoOut;


		position += block_size;
		cout<<"position+block_size: "<<position+block_size<<endl;
		cout<<"iq_data.size(): "<<iq_data.size()<<endl;
	}
}

// void stereo(const int mode,std::vector<float>& audio_data)
// {
	
// 	fmPLL(std::vector<float> PLLin, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth,std::vector<float> state, std::vector<float> &ncoOut, float feedbackI, float feedbackQ);
// }
// void frontend(const int mode,std::vector<float>& audio_data)
// {
// 	float feedbackI = 1.0;
// 	float feedbackQ = 0.0;
// 	std::vector<float> PLLin;
// 	std::vector<float> ncoOut;
// 	fmPLL(std::vector<float> PLLin, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth,std::vector<float> state, std::vector<float> &ncoOut, float feedbackI, float feedbackQ);
// }


int main()
{
	int mode = 2;
	std::vector<float> audio_data; //output audio sample vector
	
	mono(mode, audio_data);
	cout <<"size of output: "<<audio_data.size()<<endl;

	const std::string out_fname = "../data/float32samples.bin";
	writeBinData(out_fname,audio_data);

	return 0;
}