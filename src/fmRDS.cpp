#include "dy4.h"
#include "filter.h"
#include <math.h>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>

void rds_thread(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state,int mode){
  if(mode == 0 || mode == 2){
    std::vector<float> RDS_extract_ceff, extract_state;
    std::vector<float> square_coeff, square_state;
    std::vector<float> rds_lpf3k_coeff, lpf3k_state;
    std::vector<float> anti_img_coeff, anti_img_state;
    std::vector<float> rrc_coeff, rrc_state;
    std::vector<float> extract_rds,extract_rds_squared,pre_PLL_rds, post_PLL_rds, mixed,lpf_filt_rds, resample_rds, rrc_rds, rrc_coeff;
    std::vector<float> fm_demod;
    //first BPF for carrier recovery
    float RDS_freq_b = 54000;
    float RDS_freq_e = 60000; //beginning and ending frequency of rds
    float Fs = 240e3;
    //second BPF for channel extraction
    float squared_freq_b = 113500;
    float squared_freq_e = 114500;
    //PLL parameters
    float center_freq = 114000;
    float phase_adjust = 0; //this value should be obtained by plotting constalation, however not enough time to do that
    pll_state statePLL;
    statePLL.integrator = 0.0;
		statePLL.phaseEst = 0.0;
		statePLL.feedbackI = 1.0;
		statePLL.feedbackQ = 0.0;
		statePLL.trigOffset = 0.0;
		statePLL.ncoLast = 1.0;
    //3khz cutoff
    float LPF_3K_Fc = 3e3;
    //rational resampler
    int upsample;
    int downsample;
    float anti_img_cutoff;
    //clock
    std::vector<float> symbols_I;
    //mode 0 expected frequency = 23 * 2375 = 54625Hz, mode 2 expected freqency = 80750Hz, if_Fs = 240K
      if(mode == 0){
        upsample = 437;
        downsample = 1920;
        anti_img_cutoff = 54625/2;
      }
      else if(mode == 2){
        upsample = 323;
        downsample = 960;
        anti_img_cutoff = 80750/2;
      }

    //initialize vectors with size
		pre_state_extract.resize(num_taps-1);
		square_state.resize(num_taps-1);
		lpf_3k_state.resize(num_taps-1);
		anti_img_state.resize(num_taps*19-1);
		rrc_state.resize(num_taps-1);
    //Get coeeficents
		impulseResponseBPF(RDS_freq_b, RDS_freq_e, Fs, num_taps, extract_RDS_coeff);
		impulseResponseBPF(squared_freq_b, squared_freq_e, Fs, num_taps, square_coeff);
		impulseResponseLPF(Fs, LPF_3K_Fc, num_taps, rds_lpf3k_coeff);
    impulseResponseLPF(Fs*upsample, anti_img_cutoff, num_taps*upsample, anti_img_coeff);
      if(mode == 0){
        impulseResponseRRC(rrc_coeff, 54625, num_taps);
      }
      else if(mode == 2){
		    impulseResponseRRC(rrc_coeff, 80750, num_taps);
      }

      while(true)
		  {
			//Creates lock
			std::unique_lock<std::mutex> queue_lock(queue_mutex);
			//Waits until there is something in the queue
			if(data_queue.empty())
			{
				processing.wait(queue_lock);
			}
			//Assigns front of quene to vector
			fm_demod = data_queue.front();
			//Pops the value of the queue
			data_queue.pop();
			queue_lock.unlock();
			processing.notify_one();

			// RDS channel extraction
			convolutionBlock(extract_rds, fm_demod, extract_RDS_coeff, extract_state, 1.0);

			// carrier recovery
			//phase_adjust should be determined by looking at the constalation diagrams. Not enough time to look into it
			//do squaring first
      extract_rds_squared.resize(extract_rds.size(), 0.0);
      for(unsigned int i = 0; i < extract_rds.size(); i++){
        extract_rds_squared[i] = extract_rds[i] * extract_rds[i];
      }
      //now do bpf
      convolutionBlock(pre_PLL_rds, extract_rds_squared, square_coeff, square_state, 1);
      //after this do Pll
      fmPLL(post_PLL_rds, pre_PLL_rds, center_freq, 240000, 2.0, 0.0, 0.01, statePLL);

			// demodulation mixer
			//Low pass filter combined with mixer make it faster
			convolutionCombined(lpf_filt_rds, post_Pll_rds, extract_rds, rds_lpf3k_coeff, lpf3k_state, 1);

			//Resampler and pass the signal through anti imaging filter using its impluse response
			convolutionResampling(resample_rds, lpf_filt_rds, anti_img_coeff, anti_img_state, downsample, upsample);
      for(unsigned int i = 0; i<resample_rds.size(); i++){
        resamples_rds[i] *= upsample;
      }


			//RRC convolution
			convolutionBlock(rrc_rds, resample_rds, rrc_coeff, rrc_state, 1);

			//Push to thread which dedicated to the rest of data processing in RDS
			std::unique_lock<std::mutex> frame_lock(frame_mutex);
			if(data_queue.size() == 5)
			{
				processing.wait(queue_lock);
			}
			//Push onto queue
			data_queue.push(rrc_rds);
			queue_lock.unlock();
			processing.notify_one();
			//clear the data with zeros
			std::fill(extract_rds.begin(), extract_rds.end(), 0);
			std::fill(pre_Pll_rds.begin(), pre_Pll_rds.end(), 0);
			std::fill(lpf_filt_rds.begin(), lpf_filt_rds.end(), 0);
			std::fill(resample_rds.begin(), resample_rds.end(), 0);
			std::fill(rrc_rds.begin(), rrc_rds.end(), 0);
			//iterate block id
			block_id ++;

			//Will keep iterating the block till there is nothing coming in from standard in
			if(end_state && data_queue.empty())
			{
				break;
			}
		}
  }
}
//Thread which does all of the clock recovery, manchester decoding and frame sync (incomplete)
void frame_thread(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state){
	if(mode == 0 || mode == 2)
	{
		//Block id
		int block_id = 0;
		//Stuff for clock sync
		unsigned int initial_offset = 0;
		//variables for manchester decoding
		int count_0_pos =0;
		int count_1_pos = 0;
		unsigned int start_pos = 0;
		//variable for differential decoding
		int prebit = 0;
		unsigned int offset = 0;
		//initialize some vector for incoming data and state saving during processing
		std::vector<float> rrc_rds;
		std::vector<int> bit_stream;
		std::vector<int> diff_bits;
		std::vector<float> symbols_I;
		std::vector<int> potential_syndrome,block,prev_sync_bits;
		//Resize the vector
		potential_syndrome.resize(10);
		block.resize(26);
		prev_sync_bits.resize(27);
		//Poistion variables for frame sync
		unsigned int position = 0;
		unsigned int printposition = 0;
		int last_position = -1;
		//The parity matrix
		std::vector<std::vector<int>> H {{1,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},{1,0,1,1,0,1,1,1,0,0},{0,1,0,1,1,0,1,1,1,0},{0,0,1,0,1,1,0,1,1,1},{1,0,1,0,0,0,0,1,1,1},{1,1,1,0,0,1,1,1,1,1},{1,1,0,0,0,1,0,0,1,1},{1,1,0,1,0,1,0,1,0,1},{1,1,0,1,1,1,0,1,1,0},{0,1,1,0,1,1,1,0,1,1},{1,0,0,0,0,0,0,0,0,1},{1,1,1,1,0,1,1,1,0,0}, {0,1,1,1,1,0,1,1,1,0},{0,0,1,1,1,1,0,1,1,1},{1,0,1,0,1,0,0,1,1,1},{1,1,1,0,0,0,1,1,1,1}, {1,1,0,0,0,1,1,0,1,1}};
		//define syndromes
		std::vector<int> syndrome_A {1,1,1,1,0,1,1,0,0,0};
		std::vector<int> syndrome_B {1,1,1,1,0,1,0,1,0,0};
		std::vector<int> syndrome_C {1,0,0,1,0,1,1,1,0,0};
		std::vector<int> syndrome_D {1,0,0,1,0,1,1,0,0,0};

	}
}
