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
