/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, int upFactor)
{
	h.clear(); h.resize(num_taps, 0.0);

	double Norm_cutoff = Fc/(Fs/2);
	for(int i=0;i<num_taps;i++){
		if(i==((num_taps-1)/2)){
			h[i]= Norm_cutoff;
		}
		else{
			double param = PI*Norm_cutoff*(i-((float)num_taps-1.0)/2.0);
			h[i] = Norm_cutoff*sin(param)/param;
		}
		h[i]=h[i]*pow(sin(i*PI/num_taps),2)*(float)upFactor;
	}
}

void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h, int upFactor)
{
	h.clear(); h.resize(num_taps, 0.0);

	double normCenter = ((Fe+Fb)/2)/(Fs/2);
	double normPass = (Fe-Fb)/(Fs/2);

	for(int i=0;i<num_taps;i++){
		if(i==((num_taps-1)/2)){
			h[i]= normPass;
		}
		else{
			double param = PI*normPass/2*(i-((float)num_taps-1.0)/2.0);
			h[i] = normPass*sin(param)/param;
		}
		h[i]=h[i]*cos((i-(num_taps-1)/2)*PI*normCenter);
		h[i]=h[i]*pow(sin(i*PI/num_taps),2)*(float)upFactor;
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for(int n=0;n<y.size();n++){
		for(int k=0;k<h.size();k++){
			if((n-k)>=0 && (n-k)<x.size()){
				y[n]+=h[k]*x[n-k];
			}
		}
	}
}

void blockConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{
	// allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size(), 0.0);

	for(int n=0;n<y.size();n++){
		for(int k=0;k<h.size();k++){
			if((n-k)>=0){
				y[n]+=h[k]*x[n-k];
			} else {
				y[n] += h[k] * state[state.size() - (k-n)];
			}
		}
	}

	state.assign(x.end() - state.size(), x.end());
}

void fmDemodArctan(const std::vector<float> &I, const std::vector<float> &Q, float &prev_I, float &prev_Q, std::vector<float>& fm_demod) {
	fm_demod.resize(I.size());
	for(int k=0;k<I.size();k++){
		float param = (std::pow(I[k],2)+std::pow(Q[k],2));
		if(param == 0){
			fm_demod[k]=0.0;
			continue;
		}
		if(k>0){
			fm_demod[k] = ((I[k]*(Q[k]-Q[k-1]))-(Q[k]*(I[k]-I[k-1])))/param;
		}
		else{
			fm_demod[k] = (I[k]*(Q[k]-prev_Q)-Q[k]*(I[k]-prev_I))/param;
		}
	}
	prev_I = I[I.size()-1];
	prev_Q = Q[Q.size()-1];
}

void downsample(const std::vector<float> data, size_t factor, std::vector<float>& downsampled) {
    // Iterate through the data and take every 'factor' element
	downsampled.clear();
    for (int i = 0; i < data.size(); i += factor) {
        downsampled.push_back(data[i]);
    }
}

void upsample(const std::vector<float> data, size_t factor, std::vector<float> &upsampled){
    //iterate through data and insert "factor" number of zeroes in front of each data point
	upsampled.clear();
    for (int i = 0; i < data.size(); i++) {
        upsampled.push_back(data[i]);
        for(int j = factor; j > 1; j--){
            upsampled.push_back(0);
        }
    }
}

void downsampleBlockConvolveFIR(int factor, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{
    // allocate memory for the output (filtered) data
    y.clear();
    y.resize(x.size()/factor, 0.0);

    for(int n = 0; n < x.size(); n += factor){
        for(int k=0;k<h.size();k++){
            if((n-k)>=0){
                y[n/factor]+=h[k]*x[n-k];
            } else {
                y[n/factor] += h[k] * state[state.size() - (k-n)];
            }
        }
    }

    state.assign(x.end() - state.size(), x.end());
}

void resampleBlockConvolveFIR(int upFactor, int downFactor, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state)
{

	// static int debug_block = 0;

    // allocate memory for the output (filtered) data
    y.clear(); //y.resize(x.size()+h.size()-1, 0.0);
    y.resize((x.size()/(float)downFactor)*upFactor, 0.0);

	// if (debug_block < 11 && debug_block >9){
	// 	std::cout<<"xb size: "<<xb.size()<<std::endl;
	// 	std::cout<<"yb size: "<<yb.size()<<std::endl;
	// 	std::cout<<"h size: "<<h.size()<<std::endl;

	// }

    for(int n = 0; n < x.size()*(upFactor); n += downFactor){
        int phase = n % upFactor;
        for(int k = phase; k < h.size(); k += upFactor){
            if((n-k)>=0){
                y[n/downFactor]+=h[k] * x[(n-k)/upFactor];
            } else {
                y[n/downFactor] += h[k] * state[state.size() - ((k-n)/upFactor)];
            }
        }
    }

    state.assign(x.end() - state.size(), x.end());

	// debug_block++;

}
void fmPLL(const std::vector<float> &PLLin, const float freq, const float Fs, const float ncoScale, const float phaseAdjust, const float normBandwidth, std::vector<float> &ncoOut, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &trigOffset, float &nco_state){
	float Cp = 2.666;
	float Ci = 3.555;

	float Kp = normBandwidth*Cp;
	float Ki = normBandwidth*normBandwidth*Ci;

	//ncoOut.clear();
	ncoOut.resize(PLLin.size(), 0.0);

	ncoOut[0] = nco_state;
	// int trigOffset = 0;

	float errorI, errorQ, errorD;
	float trigArg;

	for (int k = 0; k < PLLin.size(); k++) {
		// phase detector
		errorI = (PLLin[k] == 0 ? 1 : PLLin[k]) * (feedbackI);
		errorQ = PLLin[k] * (-1*feedbackQ);

		if (PLLin[k] != PLLin[k]) {
		//	std::cout<<"PLLin[k] touched"<<std::endl;
		}

		// four-quadrant arctangent discriminator for phase error detection
		errorD = atan2(errorQ,errorI);

		if (errorD != errorD) {
		//	std::cout<<"errorD touched"<<std::endl;
		}

		// loop filter
		integrator += Ki*errorD;

		// update phase estimate
		phaseEst += Kp*errorD + integrator;

		// internal oscillator
		trigOffset++;
		trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;

		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		if (k == PLLin.size() - 1) {
			nco_state = cos(trigArg*ncoScale + phaseAdjust);
		} else {
			ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
		}
	}

	// for stereo only the in-phase NCO component should be returned
	// for block processing you should also return the state
	// for RDS add also the quadrature NCO component to the output
}
void delayBlock(const std::vector<float>&input_block, std::vector<float>&state_block, std::vector<float>&output_block){
	//fm_demod_size: 5120, state_block: 101-1=100
	//output_block.clear();
	
	output_block.resize(input_block.size());
	
	std::copy(state_block.begin(), state_block.end(), output_block.begin());
	std::copy(input_block.begin(), input_block.end() - state_block.size(), output_block.begin() + state_block.size());
	//state_block.assign(input_block.end() - state_block.size(), input_block.end());
	
	std::copy(input_block.end() - state_block.size(), input_block.end(), state_block.begin());

	//std::cout<<"state_block.size(): "<<state_block.size()<<std::endl;
	//std::cout<<"input_block.size(): "<<input_block.size()<<std::endl;

	//output_block.insert(output_block.end(), state_block.begin(), state_block.end());
	//output_block.insert(output_block.begin()+state_block.size(), input_block.begin(), input_block.end()-state_block.size());

	//int stateSize = state_block.size();
	//state_block.clear();
	// state_block.assign(input_block.end()-stateSize, input_block.end());
	//state_block.insert(state_block.end(), input_block.end()-stateSize, input_block.end());
}

void pointwiseMultiply(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output){
	//output.clear();
	// if(block1.size()!=block2.size()){
	// 	std::cout<<"size mismatch in multiply mixer"<<std::endl;
	// 	std::cout<<"blk1 size: "<<block1.size()<<std::endl;		
	// 	std::cout<<"blk2 size: "<<block2.size()<<std::endl;
	// }
	int size=block1.size()<block2.size()?block1.size():block2.size();
	output.resize(size);

	for(int i=0;i<size;i++){
		output[i] = block1[i]*block2[i]*2; //add gain to each element
	}
}
void pointwiseAdd(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output){
	//output.clear();
	output.resize(block1.size());
	// if(block1.size()!=block2.size()){
	// 	std::cout<<"size mismatch in add mixer"<<std::endl;
	// 	std::cout<<"blk1 size: "<<block1.size()<<std::endl;		
	// 	std::cout<<"blk2 size: "<<block2.size()<<std::endl;
	// }
	for(int i=0;i<block1.size();i++){
		output[i] = block1[i]+block2[i];
	}
}
void pointwiseSubtract(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output){
	//output.clear();
	output.resize(block1.size());
	// if(block1.size()!=block2.size()){
	// 	std::cout<<"size mismatch in subtract mixer"<<std::endl;
	// 	std::cout<<"blk1 size: "<<block1.size()<<std::endl;		
	// 	std::cout<<"blk2 size: "<<block2.size()<<std::endl;
	// }
	for(int i=0;i<block1.size();i++){
		output[i] = block1[i]-block2[i];
	}
}
void interleave(const std::vector<float>&left,const std::vector<float>&right,std::vector<float>&output){
	int size = left.size()+right.size();
	//output.clear();
	output.resize(size);
	for(int i=0;i<size;i+=2){
		output[i]=left[i/2];
	}
	for(int i=1;i<size;i+=2){
		output[i]=right[i/2];
	}
}
