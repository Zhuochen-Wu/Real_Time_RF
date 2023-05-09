#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "helper.h"
#include "dy4.h"
#include "fourier.h"
#include "iofunc.h"
#include "logfunc.h"



// separate I and Q data from incoming rf data
void iqSplit(const std::vector<float> &data_in, std::vector<float> &i_data, std::vector<float> &q_data)
{
    i_data.resize(data_in.size()/2);
    q_data.resize(data_in.size()/2);
    int j = 0;
    for(unsigned int i=0; i<data_in.size(); i+=2) {
        i_data[j] = data_in[i];
        q_data[j] = data_in[i+1];
        j++;
    }
}

void fmDemodArctan(std::vector<float> &fm_demod, std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase)
{
    float dQ;
    float dI;
    fm_demod.resize(I.size());
    for(unsigned int i =0; i<I.size(); i++){
        dQ = Q[i] - prev_phase[1];
        dI = I[i] - prev_phase[0];
        if(pow(I[i],2) + pow(Q[i],2)==0){
            fm_demod[i] = 0;
        } else {
            fm_demod[i] = (I[i]*dQ - Q[i]*dI)/(pow(I[i], (float)2) + pow(Q[i], (float)2));
        }

        prev_phase = {I[i], Q[i]};
    }
}



void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state &state){

    float Cp = 2.666;
	float Ci = 3.555;

	// gain for the proportional term
	float Kp = (normBandwidth)*Cp;
	// gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	// output array for the NCO
	ncoOut.resize(pllIn.size()+1);

    //Intiliazing variables used inside PLL to state saved values
    float integrator = state.integrator;
    float phaseEst = state.phaseEst;
    float feedbackI = state.feedbackI;
    float feedbackQ = state.feedbackQ;
    ncoOut[0] = state.ncoLast;
    float trigOffset = state.trigOffset;

    // float phaseAdjust = 0.0;

	for (unsigned int k=0; k<pllIn.size(); k++)
    {
		float errorI = pllIn[k] * (+feedbackI);  // complex conjugate of the
		float errorQ = pllIn[k] * (-feedbackQ); // feedback complex exponential

    	// four-quadrant arctangent discriminator for phase error detection
		float errorD = atan2(errorQ, errorI);

		// loop filter
		integrator += Ki*errorD;

		// update phase estimate
		phaseEst += Kp*errorD + integrator;

		// internal oscillator
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }

    //Update State Variables so they are saved to the struct object in main
    state.integrator = integrator;
    state.phaseEst = phaseEst;
    state.feedbackI = feedbackI;
    state.feedbackQ = feedbackQ;
    state.ncoLast= ncoOut[ncoOut.size()-1];
    state.trigOffset= trigOffset + pllIn.size();

    //Resize to return 1:end of array
    ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);

}
