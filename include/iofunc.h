/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void printRealVector(const std::vector<float> &);

void printComplexVector(const std::vector<std::complex<float>> &);

void readBinData(const std::string, std::vector<float> &);

void writeBinData(const std::string, const std::vector<float> &);

void readRawData(const std::string in_fname, std::vector<uint8_t> &raw_data);

void convertRaw(const std::vector<uint8_t> raw_data, std::vector<float> &iq_data);

void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right);

void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data);

#endif // DY4_IOFUNC_H
