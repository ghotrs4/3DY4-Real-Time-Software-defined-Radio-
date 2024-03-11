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
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void blockConvolveFIR(std::vector<float> &y, const std::vector<float> x, const std::vector<float> h, std::vector<float> &state, int position, int block_size);
void fmDemodArctan(std::vector<float> I, std::vector<float> Q, float &prev_I, float &prev_Q, std::vector<float>& fm_demod);
void downsample(const std::vector<float> data, size_t factor, std::vector<float>& downsampled);
void upsample(const std::vector<float> data, size_t factor, std::vector<float> &upsampled);

#endif // DY4_FILTER_H