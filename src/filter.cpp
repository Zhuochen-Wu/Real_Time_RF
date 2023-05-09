/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float norm_cutoff = Fc/(Fs/2);
	float Ndiv = (num_taps-1)/2;
	for(int i=0; i<num_taps; i++){
		if (i == Ndiv){
			h[i] = norm_cutoff;
		} else {
			h[i] = norm_cutoff*sin(PI*norm_cutoff*(i-Ndiv))/(PI*norm_cutoff*(i-Ndiv));
		}
		h[i] =h[i]*(float)pow(sin(i*PI/num_taps),2.0);
	}
}


// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
//this function is not utilized in our code as it is too slow
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for (unsigned int n=0;n<y.size();n++){
  	for (unsigned int k=0;k<h.size();k++){
   		if ((n-k)>=0 and n-k < x.size()){
    		y[n]=y[n]+h[k]*x[n-k];
			}
		}
	}
}
//do convolution with resampling at the same time
void convolutionResampling(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int decim, int upsample)
{
	y.resize(x.size()*upsample/decim);

	int s = state.size()-1;
	int count;

	for(unsigned int i=0; i < y.size(); i++){

		count = 0;

		// iterate through h, the step size is up sampling factor
		for(unsigned int j=0; j < h.size(); j += upsample) {
			if(decim*i-j >= 0 && decim*i-j < x.size()*upsample){
				y[i] += x[(decim*i-j)/upsample]*h[j];
			}
			else
			{
				y[i] += state[s-count] * h[j];
				count++;
			}
		}
	}

	auto first = (x.size()-(h.size()-1));
	auto last = x.end();

	state = std::vector<float>(x.begin()+first, last);
}





void impulseResponseRRC(std::vector<float> &impulseResponseRRC, float Fs, float N_taps){
  /*Root raised cosine (RRC) filter

	Fs  		sampling rate at the output of the resampler in the RDS path
				sampling rate must be an integer multipler of 2375
				this integer multiple is the number of samples per symbol

	N_taps  	number of filter taps
  */
  float t;
  //duration for each symbol - do NOT be changed for RDS!
  float T_symbol = 1/2375.0;
  //roll-off factor (greater than 0 and smaller than 1)
  float beta = 0.90;
  //the RRC inpulse response that will be computed in this function
  impulseResponseRRC.clear();
	impulseResponseRRC.resize(N_taps);
  for (int k=0; k<N_taps; k++){
    //we ignore the 1/T_symbol scale factor
    t = ((k-N_taps)/2.0)/Fs;
    if (t==0.0){
      impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
    }
    else if ((t==(-T_symbol/(4*beta))) || (t==T_symbol/(4*beta))){
      impulseResponseRRC[k] = (beta/sqrt(2))*(((1+2/PI)*(sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
    }
    else {impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/ (PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);

    }
  }
}

// filter both I and Q also do downsampling in rf front end
void convolutionBlockIQ(std::vector<float> &I, std::vector<float> &Q, const std::vector<float> &I_data, const std::vector<float> &Q_data, const std::vector<float> &h,std::vector<float> &stateI, std::vector<float> &stateQ,int decim){
	I.resize(I_data.size()/decim);
	Q.resize(Q_data.size()/decim);
	int d = 0;
	int s = stateI.size()-1;
	int count;
	 for(unsigned int n = 0 ; n < I_data.size(); n+=decim){
		 count=0;
		 for(unsigned int k = 0 ; k < h.size(); k++){
			if (n-k >= 0 && n-k <I_data.size()){
				I[d] += h[k]*I_data[n-k];
				Q[d] += h[k]*Q_data[n-k];
			}else{
				I[d] += h[k]*stateI[s-count];
				Q[d] += h[k]*stateQ[s-count];
				count++;
			}
		}
		d++;
	}

	auto first = (I_data.size()-(h.size()-1));

	stateI = std::vector<float>(I_data.begin()+first, I_data.end());
	stateQ = std::vector<float>(Q_data.begin()+first, Q_data.end());

}
//do convolution with down sample only
void convolutionBlock(std::vector<float> &y,const std::vector<float> &x, const std::vector<float> &h,std::vector<float> &state,const int decim){
	y.resize(x.size()/decim);
	int y_original = x.size();
	int d = 0;
	int s = state.size()-1;
	int count;
	 for(int n = 0 ; n < y_original; n+=decim){
		 count=0;
		 for(unsigned int k = 0 ; k < h.size(); k++){
			if (n-k >= 0 && n-k <x.size()){
				y[d] += h[k]*x[n-k];
			}else{
				y[d] += h[k]*state[s-count];
				count++;
			}
		}
		d++;
	}

	auto first = (x.size()-(h.size()-1));
	auto last = x.end();

	state = std::vector<float>(x.begin()+first, last);
}
//band pass filter for stereo carrier recovery and channel extraction
void impulseResponseBPF(float Fb, float Fe, float Fs, int n_taps, std::vector<float> &h){

	h.resize(n_taps, 0.0);
	float norm_center = ((Fe+Fb)/2)/(Fs/2);
        float norm_pass = (Fe-Fb)/(Fs/2);
	float n_half = (n_taps-1)/2;
	for(int i=0; i<n_taps; i++){
		if(i == n_half){
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass * sin(PI*(norm_pass/2)*(i-n_half)) / (PI*(norm_pass/2)*(i-n_half));
		}
		h[i] = h[i] * cos(i*PI*norm_center);
        h[i] = h[i] * pow(sin((i*PI)/n_taps), 2);
	}
}
//combine convolution, downsampling and mixer
void convolutionCombined(std::vector<float> &y, const std::vector<float> &x,const std::vector<float> &x1 , const std::vector<float> &h, std::vector<float> &state, const int &decim_num)
{
	//Creates vector for down sampled data
	y.resize(x.size()/decim_num);
	//Loops through vlaues to do convoloution
  int count;
  for (unsigned int n = 0; n < y.size(); n++)
	{
      count = 0;
      for (unsigned k = 0; k < h.size(); k++)
		  {
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
			{
				//do down sampling and mixing at the same time
        y[n] += x[decim_num*n-k]*x1[decim_num*n-k]*h[k]*2;
			}
			//Previous state data
			else
			{
        y[n] += state[state.size()-1-count]*h[k];
				count += 1;
			}
    }
  }
	//Assigns next state value
	for(unsigned int i = 0; i < state.size(); i++){
		state[i] = x[x.size()-state.size()-1+i]*x1[x.size()-state.size()-1+i];
	}
}
