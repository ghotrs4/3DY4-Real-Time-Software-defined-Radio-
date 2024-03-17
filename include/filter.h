/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &, int);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void blockConvolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int position, int block_size);
void fmDemodArctan(std::vector<float> &I, std::vector<float> &Q, float &prev_I, float &prev_Q, std::vector<float>& fm_demod);

void downsample(const std::vector<float> data, size_t factor, std::vector<float>& downsampled);
void upsample(const std::vector<float> data, size_t factor, std::vector<float> &upsampled);

void downsampleBlockConvolveFIR(int factor, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int position, int block_size);
void resampleBlockConvolveFIR(int upFactor, int downFactor, std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int position, int block_size);
void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h, int upFactor);
void fmPLL(const std::vector<float> &PLLin, const float freq, const float Fs, const float ncoScale, const float phaseAdjust, const float normBandwidth, std::vector<float> &ncoOut, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst);
void pointwiseMultiply(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output);
void delayBlock(const std::vector<float>&input_block, std::vector<float>&state_block, std::vector<float>&output_block);
void pointwiseAdd(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output);
void pointwiseSubtract(const std::vector<float>&block1,const std::vector<float>&block2,std::vector<float>&output);
void interleave(const std::vector<float>&left,const std::vector<float>&right,std::vector<float>&output);

#endif // DY4_FILTER_H