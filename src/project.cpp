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

#include <chrono>

using namespace std;

	using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    
    

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

	switch(mode) {
		case 0: //output Fs = 48k
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 24e3;
			break;
		case 1://output Fs = 36k
			rf_Fs = 1.44e6;
			audio_Fs = 288e3;
			rf_decim = 5;
			audio_decim = 8;
			audio_upsample = 1;
			audio_Fc = 18e3;
			break;
		case 2://output Fs = 44.1k
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 800;
			audio_upsample = 147;
			audio_taps *= audio_upsample;
			audio_Fs *= audio_upsample;
			audio_Fc = 22.05e3;
			break;
		case 3://output Fs = 44.1k
			rf_Fs = 1.92e6;
			audio_Fs = 128e3;
			rf_decim = 15;
			audio_decim = 1280;
			audio_upsample = 441;
			audio_taps *= audio_upsample;
			audio_Fs *= audio_upsample;
			audio_Fc = 22.05e3;
			break;
		default:
			rf_Fs = 2.4e6;
			audio_Fs = 240e3;
			rf_decim = 10;
			audio_decim = 5;
			audio_upsample = 1;
			audio_Fc = 24e3;
	}

	const std::string in_fname = "../data/1920.raw";
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

	std::vector<float> i_block;
	std::vector<float> q_block;

	std::vector<float> i_downsampled;
	std::vector<float> q_downsampled;

	std::vector<float> audio_filt;
	std::vector<float> audio_block;
	cout<<"block size"<<block_size<<endl;
	cout<<"begin debug: ..."<<endl;
	cout<<"rf_coeff sample 0: "<<rf_coeff[0]<<endl;
	
	
    
    /*
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    
    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    
    auto t5 = high_resolution_clock::now();
    auto t6 = high_resolution_clock::now();
    
    auto t7 = high_resolution_clock::now();
    auto t8 = high_resolution_clock::now();
    
    auto tbeg = high_resolution_clock::now();
    
    auto tend = high_resolution_clock::now();
    */
    
	while (position+block_size<iq_data.size()) {//touched
		//tbeg = high_resolution_clock::now();
		
		cout<<"block number: "<<position/block_size<<endl;

		// for(int i=0;i<5;i++){
		// 	cout<<"i block samples: "<<i_samples[i+position/2]<<endl;
		// }
		// for(int i=0;i<5;i++){
		// 	cout<<"q block samples: "<<q_samples[i+position/2]<<endl;
		// }
		
		//t1 = high_resolution_clock::now();
		
		//blockConvolveFIR(i_block, i_samples, rf_coeff, i_state_rf, position/2, block_size/2);
		//blockConvolveFIR(q_block, q_samples, rf_coeff, q_state_rf, position/2, block_size/2);
		
		downsampleBlockConvolveFIR(rf_decim, i_downsampled, i_samples, rf_coeff, i_state_rf, position/2, block_size/2);
		downsampleBlockConvolveFIR(rf_decim, q_downsampled, q_samples, rf_coeff, q_state_rf, position/2, block_size/2);
		

		
		//t2 = high_resolution_clock::now();
		
		
		// for(int i=0;i<5;i++){
		// 	cout<<"i block convolved: "<<i_block[i]<<endl;
		// }
		// for(int i=0;i<5;i++){
		// 	cout<<"q block convolved: "<<q_block[i]<<endl;
		// }
		
		//t3 = high_resolution_clock::now();
		//downsample(i_block, rf_decim, i_downsampled);
		//downsample(q_block, rf_decim, q_downsampled);
		
		//t4 = high_resolution_clock::now();
		
		
		// for(int i=0;i<5;i++){
		// 	cout<<"i down samples: "<<i_downsampled[i]<<endl;
		// }
		// for(int i=0;i<5;i++){
		// 	cout<<"q down samples: "<<q_downsampled[i]<<endl;
		// }
		
		//t5 = high_resolution_clock::now();
		fmDemodArctan(i_downsampled, q_downsampled, prev_I, prev_Q, fm_demod);
		//t6 = high_resolution_clock::now();
		
		// for(int i=0;i<5;i++){
		// 	cout<<"fmdemod samples: "<<fm_demod[i]<<endl;
		// }
		cout<<"size of fmdemod: "<<fm_demod.size()<<endl;
		// cout<<"size of block size: "<<block_size/rf_decim<<endl;
		// cout << "wtf" << endl;
		// for(int i=0;i<5;i++){
		// 	cout<<"audio filt samples: "<<audio_filt[i]<<endl;
		// }

		// for(int i=0;i<5;i++){
		// 	cout<<"audio block samples: "<<audio_block[i]<<endl;
		// }
		// downsampleBlockConvolveFIR(audio_decim, audio_block, fm_demod, audio_coeff, state_audio, 0, fm_demod.size());
		
		//t7 = high_resolution_clock::now();
		resampleBlockConvolveFIR(audio_upsample, audio_decim, audio_block, fm_demod, audio_coeff, state_audio, 0, fm_demod.size());
		//t8 = high_resolution_clock::now();
		
		if (position != 0) {
			audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());

		}

		position += block_size;
		cout<<"position+block_size: "<<position+block_size<<endl;
		cout<<"iq_data.size(): "<<iq_data.size()<<endl;
	
//		tend = high_resolution_clock::now();
		
	}
	
	
	/* Getting number of milliseconds as a double. */
	//duration<double, std::milli> ms_double1 = t2 - t1;
	
	//duration<double, std::milli> ms_double2 = t4 - t3;
	
	
	//duration<double, std::milli> ms_double3 = t6 - t5;
	
	
	//duration<double, std::milli> ms_double4 = t8- t7;
	
	//duration<double, std::milli> ms_double5 = tend- tbeg;
	
	
	//Prints times
	//std::cout << "iq resampleblockConvoleFir time = " << ms_double1.count() << "ms\n";
	
	
	//std::cout << "iq downsampling time = " << ms_double2.count() << "ms\n";
	
	//std::cout << "fmdemod Arctan time = " << ms_double3.count() << "ms\n";
	//std::cout << "resample time = " << ms_double4.count() << "ms\n";
	
	cout<<"block size: " << block_size<<endl;
	
	
}

int main()
{
	auto timeStart = high_resolution_clock::now();
	int mode = 3;
	std::vector<float> audio_data; //output audio sample vector
    
    
    
	mono(mode, audio_data);
	

	
	//std::cout << "total time = " << ms_double5.count() << "ms\n";
	
	cout <<"size of output: "<<audio_data.size()<<endl;

	const std::string out_fname = "../data/float32samples.bin";
	writeBinData(out_fname,audio_data);
	
	auto timeEnd = high_resolution_clock::now();
	
	duration<double, std::milli> ms_double6 = timeEnd- timeStart;
	
	std::cout << "total time = " << ms_double6.count() << "ms\n";

	return 0;
}
