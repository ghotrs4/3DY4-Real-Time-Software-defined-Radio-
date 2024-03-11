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
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	h.clear(); h.resize(num_taps, 0.0);

	double Norm_cutoff = Fc/(Fs/2);
	for(int i=0;i<num_taps;i++){
		if(i==((num_taps-1)/2)){
			h[i]= Norm_cutoff;
		}
		else{
			double param = PI*Norm_cutoff*(i-(num_taps-1)/2);
			h[i] = Norm_cutoff*sin(param)/param;
		}
		h[i]=h[i]*pow(sin(i*PI/num_taps),2);
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

void blockConvolveFIR(std::vector<float> &y, const std::vector<float> x, const std::vector<float> h, std::vector<float> &state, int position, int block_size)
{
	// allocate memory for the output (filtered) data
	y.clear(); //y.resize(x.size()+h.size()-1, 0.0);

	std::vector<float> xb;
	std::vector<float> yb;


	xb = std::vector<float>(x.begin() + position, x.begin() + position + block_size); // new block
	yb.clear(); // clear output
	yb.resize(xb.size(), 0.0);

	for(int n=0;n<yb.size();n++){
		for(int k=0;k<h.size();k++){
			if((n-k)>=0){
				yb[n]+=h[k]*xb[n-k];
			} else {
				yb[n] += h[k] * state[state.size() - (k-n)];
			}
		}
	}

	if (state.size() > block_size) {
		state.insert(state.begin(), state.begin() + block_size, state.end()); // left shift to make room
		state.insert(state.end() - block_size, xb.begin(), xb.end()); // put whole block into state
	} else {
		state = std::vector<float>(xb.end() - state.size(), xb.end());
	}

	y=yb;
}

void fmDemodArctan(std::vector<float> I, std::vector<float> Q, float &prev_I, float &prev_Q, std::vector<float>& fm_demod) {
	fm_demod.resize(I.size());
	for(int k=0;k<I.size();k++){
		if(k>0){
			fm_demod[k] = (I[k]*Q[k-1]-Q[k]*I[k-1])/(pow(I[k],2)+pow(Q[k],2));
		}
		else{
			fm_demod[k] = (I[k]*prev_Q-Q[k]*prev_I)/(pow(I[k],2)+pow(Q[k],2));
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
void downsampleBlockConvolveFIR(size_t factor, std::vector<float> &y, const std::vector<float> x, const std::vector<float> h, std::vector<float> &state, int position, int block_size)
{
    // allocate memory for the output (filtered) data
    y.clear(); //y.resize(x.size()+h.size()-1, 0.0);

    std::vector<float> xb;
    std::vector<float> yb;


    xb = std::vector<float>(x.begin() + position, x.begin() + position + block_size); // new block
    yb.clear(); // clear output
    yb.resize(xb.size()/factor, 0.0);
    
    int g = 0;

    for(int n = 0; n < xb.size(); n += factor){
        for(int k=0;k<h.size();k++){
            if((n-k)>=0){
                yb[g]+=h[k]*xb[n-k];
            } else {
                yb[g] += h[k] * state[state.size() - (k-n)];
            }
        }
        g++;
    }

    if (state.size() > block_size) {
        state.insert(state.begin(), state.begin() + block_size, state.end()); // left shift to make room
        state.insert(state.end() - block_size, xb.begin(), xb.end()); // put whole block into state
    } else {
        state = std::vector<float>(xb.end() - state.size(), xb.end());
    }

    y=yb;
}