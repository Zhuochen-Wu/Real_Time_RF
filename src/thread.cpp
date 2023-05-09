#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>
#include "helper.h"
#include <functional>
#include <thread>
#include <condition_variable>
#include <queue>
#include <mutex>


void data(pll_state &samples){
	samples.integrator = 0.0;
	samples.phaseEst = 0.0;
	samples.feedbackI = 1.0;
	samples.feedbackQ = 0.0;
	samples.trigOffset = 0.0;
	samples.ncoLast = 1.0;
}// initial state of PLL

void front_end_producer(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state,int mode){
	int rf_Fs = 2400e3;
	int rf_Fc = 100000;
	int rf_taps = 151;
	int rf_decim = 10;
	unsigned int block_size = 102400;

	if(mode == 0){
		rf_Fs = 2400e3;
		rf_decim = 10;
	    block_size = 102400; //block size = 1024 * rf_decim * audio_decim * 2
	} else if(mode == 1){
		rf_Fs = 2304e3;
		rf_decim = 8;
		block_size = 98304;
	}else if(mode == 2){
		rf_Fs = 2400e3;
		rf_decim = 10;
		block_size = 3010560; //which is 1024 * 147 * 10 * 2
	}else if(mode == 3){
		rf_Fs = 960e3;
		rf_decim = 3;
		block_size = 2709504; //which is 1024 * 441 * 3200 / 100
	}

	unsigned int block_id = 0;
	std::vector<float> iq_data, i_data, q_data,iq_ds_coeff;
	std::vector<float> i_ds, q_ds,i_state,q_state;


	std::vector<float> prev_phase, fm_demod;
	std::vector<float> audio_coeff;


	//Sets some inital values
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);

	prev_phase.resize(2,0.0);
	i_state.resize(rf_taps-1,0.0);
	q_state.resize(rf_taps-1,0.0);

	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, iq_ds_coeff);

	while(true)
	{

		//read data from std in
		readStdInBlock(block_size, block_id, iq_data);
		//separate i and q data from rf
		iqSplit(iq_data, i_data, q_data);
		//filter the i q data with down sampling
		convolutionBlockIQ(i_ds, q_ds, i_data, q_data, iq_ds_coeff, i_state, q_state,rf_decim);
		//get fm demodulated data from down sampled i and q data
		fmDemodArctan(fm_demod,i_ds,q_ds,prev_phase);
		//load the queue with data from front end
		std::unique_lock<std::mutex> queue_lock(queue_mutex);
		if(data_queue.size() == 5){
			processing.wait(queue_lock);
		}
		data_queue.push(fm_demod);
		block_id++;
		queue_lock.unlock();
		processing.notify_one();
        //clear the vector for i and q down sampled for next processing
		std::fill(i_ds.begin(), i_ds.end(), 0);
		std::fill(q_ds.begin(), q_ds.end(), 0);
        //continuesly running until there is no data coming from std in
		if((std::cin.rdstate()) != 0)
		{
			end_state = true;
			break;
		}
	}
}

void stereo_consumer(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state,int mode){
	int rf_taps = 151;
	float upSampling_factor = 1;
	int if_Fs = 240000;
	int audio_Fc = 16000;
	int audio_decim = 5;
	int pilot_b = 18.5e3;
	int pilot_e = 19.5e3;
	int stereo_b = 22e3;
	int stereo_e = 54e3;
	int num_taps = 151;
	int pilot_frequency = 19e3;
	unsigned int block_size = 102400;

	if(mode == 0){
		audio_decim = 5;
		if_Fs = 240000;
		block_size = 102400;
	} else if(mode == 1){
		audio_decim = 6;
		if_Fs = 288000;
		block_size = 98304;
	}else if(mode == 2){
		audio_decim = 800;
		if_Fs = 35280000; //which is 240000 * 147
		upSampling_factor = 147;
		num_taps = 151 * 147;
		block_size = 3010560; //which is 1024 * 147
	}else if(mode == 3){
		audio_decim = 3200;
		if_Fs = 141120000; //which is 320000 * 441
		upSampling_factor = 441;
		block_size = 2709504; //which is 1024 * 441 * 3 * 2
		num_taps = 200 * 441; //increase the number of taps from 151 * 441 to 200 * 441, improve the output quality (noise reduced and volume is much higher)
	}
    //declare vectors
	std::vector<float> iq_data, i_data, q_data,iq_ds_coeff,i_state, q_state;
	std::vector<float> i_ds, q_ds;

    std::vector<float> state_recovery, state_extraction, stereo_state;
    std::vector<float> pilot_coeff,bpf_recovery,recovery_coeff,stereo_recovered;
    std::vector<float> recovery_pll;

	std::vector<float> prev_phase,fm_demod;
	std::vector<float> audio_coeff, audio_state,audio_block, audio_filter;
	std::vector<short int> audio_data;
	std::vector<float> stereo_filt;



    pll_state statePLL;
    data(statePLL);

    std::vector<float> mixed;
		if(mode == 0){
			mixed.resize(5120,0.0);
		}
		else if (mode == 1){
			mixed.resize(6144, 0.0); //block size / (2* 8)
		}
	  else if(mode == 2){
		mixed.resize(block_size/10, 0.0);
		}
		else if(mode == 3){
			mixed.resize(block_size/3, 0.0);
		}



	//initialize the size of vectors
	prev_phase.resize(2,0.0);
	audio_state.resize(num_taps-1,0.0);
	audio_state.resize(num_taps-1,0.0);
    state_recovery.resize(rf_taps-1,0.0);
    state_extraction.resize(rf_taps-1,0.0);
    stereo_state.resize(rf_taps-1,0.0);
    std::vector<float> left,right;
    std::vector<float> audio_comb;

	// LPF to extract audio coeff
	impulseResponseLPF(if_Fs, audio_Fc, num_taps, audio_coeff);
	// BAND PASS COEFF FOR PILOT AND STEREO RECOVERY
	impulseResponseBPF(pilot_b, pilot_e,if_Fs,151,pilot_coeff);
	impulseResponseBPF(stereo_b, stereo_e,if_Fs,151,recovery_coeff);


	while(true)
	{

		// UNLOCK WHEN DATA ON QUEUE
		std::unique_lock<std::mutex> queue_lock(queue_mutex);
		if(data_queue.empty()) {
			processing.wait(queue_lock);
		}
		fm_demod = data_queue.front();
		data_queue.pop();

		queue_lock.unlock();
		processing.notify_one();

		convolutionBlock(bpf_recovery, fm_demod, pilot_coeff, state_recovery, 1);
		fmPLL(recovery_pll, bpf_recovery, pilot_frequency, if_Fs/upSampling_factor, 2.0, 0.0, 0.01, statePLL);
		convolutionBlock(stereo_recovered, fm_demod, recovery_coeff, state_extraction, 1);


		if (mode <= 1) {
			convolutionBlock(audio_block, fm_demod, audio_coeff, audio_state, audio_decim); //get mono first
			for (unsigned int i = 0; i <stereo_recovered.size();i++){
            	mixed[i] = stereo_recovered[i]*recovery_pll[i]*2;
        }
			convolutionBlock(stereo_filt, mixed, audio_coeff, stereo_state, audio_decim);
		} else {
			convolutionResampling(audio_block,fm_demod,audio_coeff, audio_state,audio_decim,upSampling_factor);
			for (unsigned int i = 0; i <stereo_recovered.size();i++){
            	mixed[i] = stereo_recovered[i]*recovery_pll[i];
        }
			convolutionResampling(stereo_filt,mixed,audio_coeff,stereo_state,audio_decim,upSampling_factor);
		}

		// INTERLEAVING LEFT AND RIGHT CHANNELS
		audio_comb.resize(2*audio_block.size());
		for(unsigned int i = 0 ; i<audio_block.size(); i++){
			audio_comb[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
			audio_comb[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
		}

	    //end of combiner, now write to std out
		audio_data.resize(audio_comb.size());
		for(unsigned int l=0; l<audio_comb.size(); l++)
		{
			if(std::isnan(audio_comb[l])) {
				audio_data[l] =0;
			} else {
				audio_data[l] = static_cast<short int>(audio_comb[l] *16384*upSampling_factor);
			}
		}

		//WRITE DATA OUT
		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		//FILL WITH 0s TO CLEAR OLD DATA
		std::fill(audio_block.begin(), audio_block.end(), 0);
        std::fill(stereo_filt.begin(), stereo_filt.end(), 0);
        mixed.clear();

		//Will keep iterating the block till there is nothing coming in from producer
		if(end_state && data_queue.empty())
		{
			break;
		}
	}
}
