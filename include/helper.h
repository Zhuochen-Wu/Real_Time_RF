#ifndef DY4_HELPER_H
#define DY4_HELPER_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

struct pll_state{
		float integrator, phaseEst, feedbackI, feedbackQ, trigOffset,ncoLast;
	};

void iqSplit(const std::vector<float> &data_in, std::vector<float> &i, std::vector<float> &q);
void fmDemodArctan(std::vector<float> &fm_demod, std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase);
void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, pll_state &state);

#endif
