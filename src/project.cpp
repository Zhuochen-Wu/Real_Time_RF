/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

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
#include "thread.h"

int modeSelection(int argc, char *argv[]){
	int mode = 0;
	if (argc < 2){
		std::cerr << "no mode selected, will run in mode 0 by default..." << std::ends;
	}
	else if (argc == 2){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "invalid mode number detected, terminating..." << std::ends;
			exit(1);
		}
	}
	else{
		std::cerr << "please read document to properly select mode" << std::ends;
	}
	return mode;
}

int main(int argc, char* argv[])
{
	int mode = modeSelection(argc, argv);

	std::queue<std::vector<float>> data_queue;
	std::mutex queue_mutex;
	std::condition_variable processing;
	bool end_state = false;
	std::thread producer = std::thread(front_end_producer, std::ref(data_queue), std::ref(queue_mutex), std::ref(processing), std::ref(end_state),mode);
	std::thread consumer= std::thread(stereo_consumer, std::ref(data_queue), std::ref(queue_mutex), std::ref(processing), std::ref(end_state),mode);

	producer.join();
	consumer.join();


	return 0;
}
