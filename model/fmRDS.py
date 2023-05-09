#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
from fmPll import fmPll 
from fmRRC import impulseResponseRootRaisedCosine

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 240e3
audio_decim = 5

rds_b = 54e3
rds_e = 60e3
# add other settings for audio, like filter taps, ...

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is normalized between -1 and +1 and interleaved
	in_fname = "../data/my_samples_u8.raw"
	iq_data = np.fromfile(in_fname, dtype='uint8')
	iq_data = (iq_data -128.0)/128.0
	block_size = 1024 * rf_decim * audio_decim * 2
	# SHORTEN SAMPLES FOR QUICKER RUN TIME
	iq_data=iq_data[:15*307200]
	print("Read raw RF data from \"" + in_fname + "\" in float32 format")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, \
							rf_Fc/(rf_Fs/2), \
							window=('hann'))

	# coefficients for the filter to extract mono audio
	audio_coeff = np.array([]) # to be updated by you during in-lab

	# set up drawing
	fig, (ax0, ax1, ax2,ax3) = plt.subplots(nrows=4)
	fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is in KB and
	# a multiple of decimation factors
	# block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_phase = 0

    # RDS Extraction Coefficients
	rds_coeff = signal.firwin(rf_taps, [rds_b/(audio_Fs/2),rds_e/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")

	# Square bandpass

	square_coeff = signal.firwin(rf_taps, [113.5e3/(audio_Fs/2),114.5e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")

	#print('Processing block ' + str(block_count))

	i_filt= signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
	q_filt= signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

	# downsample the FM channel
	i_ds = i_filt[::rf_decim]
	q_ds = q_filt[::rf_decim]
	# FM demodulator

	fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

	# EXTRACTING RDS CHANNEL
	rds_data = signal.lfilter(rds_coeff, 1.0, fm_demod)

	ax0.psd(rds_data, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax0.set_ylabel('PSD (db/Hz)')
	ax0.set_title('Extracted RDS')

	#SQUARE RDS DATA
	squared_data = np.square(rds_data)

	#BANDPASS FILTER USING SQUARED COEFF FROM 113.5k to 114.5k 
	square_filtered = signal.lfilter(square_coeff, 1.0, squared_data)
	state_Pll =[0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
	phase_adjust = -1.14
	# PLL 
	recovery_I, recovery_Q, state_Pll = fmPll(square_filtered,114e3, audio_Fs, state_Pll, 0.5, phase_adjust)

	ax1.psd(recovery_I, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax1.set_ylabel('PSD (db/Hz)')
	ax1.set_title('Recovery pll')

	# MIXING WITH EXTRACTED RDS
	mixed_rds = np.multiply(rds_data, recovery_I[0:len(rds_data)])*2
	mixed_rds_Q = np.multiply(rds_data, recovery_Q[0:len(rds_data)])*2


	##LPF TO CUTOFF THE HIGHBANDS 3k

	lpf_rds = signal.firwin(rf_taps,3000/(240000/2), window=('hann'))
	lpf_filt = signal.lfilter(lpf_rds,1.0,mixed_rds)
	lpf_filtQ = signal.lfilter(lpf_rds,1.0,mixed_rds_Q)

	ax2.psd(lpf_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax2.set_ylabel('PSD (db/Hz)')
	ax2.set_title('Mixed rds Inphase')
	plt.show()
	
	#do upsampling appears to be slow since this is not optimized wit algorithm
	up = 437
	down = 1920
	up_rds = np.zeros(len(lpf_filt)*up)
	up_rdsQ = np.zeros(len(lpf_filtQ)*up)

	for i in range(len(lpf_filt)):
		up_rds[i*up] = lpf_filt[i]
		up_rdsQ[i*up] = lpf_filtQ[i]


	post_up_coeff = signal.firwin(rf_taps, (57000/2)/((240000*up)/2), window=('hann'))

	post_up = signal.lfilter(post_up_coeff,1.0,up_rds)
	post_up_q = signal.lfilter(post_up_coeff,1.0,up_rdsQ)


	resampled =  post_up[::down]*up
	resampled_q =  post_up_q[::down]*up

	##do RRC

	rrc_coeff = impulseResponseRootRaisedCosine(57000, 151)

	rrc_filt = signal.lfilter(rrc_coeff,1.0,resampled)
	rrc_filtQ = signal.lfilter(rrc_coeff,1.0,resampled_q)

	# TODO: PROPER ALGORITHM FOR RATIONAL RESAMPLER UNCOMPLETE, HARDCODED AS A TEST 
	offset = 10

	symbols_I = rrc_filt[offset::24]
	symbols_Q = rrc_filtQ[offset::24]


	fig, (phase) = plt.subplots(nrows=1)
	fig.subplots_adjust(hspace = 1.0)
	phase.scatter(symbols_I, symbols_Q, s=10)
	phase.set_ylim(-1.25, 1.25)
	plt.show()
	
	# MANCHESTER DECODING
	out_stream = np.zeros(len(symbols_I)//2)
	j = 0
	print(len(symbols_I))
	for i in range(0,len(symbols_I),2):
		if(symbols_I[i]>symbols_I[i+1]):
			out_stream[j] = 1
		j+=1

	# DIFFERENTIAL DECODING

	diff_stream = np.zeros(len(out_stream)-1)
	prev = out_stream[0]

	for i in range(len(out_stream)-1):
		diff_stream[i] = (prev and not out_stream[i+1]) or (not prev and out_stream[i+1])
		prev = out_stream[i+1]

	# syndrome and frame synchronization 
	parity_matrix = np.matrix([
							  [1,0,0,0,0,0,0,0,0,0],
							  [0,1,0,0,0,0,0,0,0,0],
							  [0,0,1,0,0,0,0,0,0,0],
							  [0,0,0,1,0,0,0,0,0,0],
							  [0,0,0,0,1,0,0,0,0,0],
							  [0,0,0,0,0,1,0,0,0,0],
							  [0,0,0,0,0,0,1,0,0,0],
							  [0,0,0,0,0,0,0,1,0,0],
							  [0,0,0,0,0,0,0,0,1,0],
							  [0,0,0,0,0,0,0,0,0,1],
							  [1,0,1,1,0,1,1,1,0,0],
							  [0,1,0,1,1,0,1,1,1,0],
							  [0,0,1,0,1,1,0,1,1,1],
							  [1,0,1,0,0,0,0,1,1,1],
							  [1,1,1,0,0,1,1,1,1,1],
							  [1,1,0,0,0,1,0,0,1,1],
							  [1,1,0,1,0,1,0,1,0,1],
							  [1,1,0,1,1,1,0,1,1,0],
							  [0,1,1,0,1,1,1,0,1,1],
							  [1,0,0,0,0,0,0,0,0,1],
							  [1,1,1,1,0,1,1,1,0,0],
							  [0,1,1,1,1,0,1,1,1,0],
							  [0,0,1,1,1,1,0,1,1,1],
							  [1,0,1,0,1,0,0,1,1,1],
							  [1,1,1,0,0,0,1,1,1,1],
							  [1,1,0,0,0,1,1,0,1,1]
							  ])

	# syndrome_matrix = np.matrix([
	# 							[1,1,1,1,0,1,1,0,0,0], ## A
	# 							[1,1,1,1,0,1,0,1,0,0], ## B
	# 							[1,0,0,1,0,1,1,1,0,0], ## C
	# 							[1,1,1,1,0,0,1,1,0,0], ## C prime
	# 							[1,0,0,1,0,1,1,0,0,0]  ## D
	# 							])

	## 1x10 matrix
	synced = 0
	pos = 0
	syndrome = ""
	while synced == False:
		## take 26 bits of bitstream, slide down by 1 if not yet synced
		check_block = diff_stream[pos:pos+26]
		result = np.zeros(10)
		## create 1x10 matrix via bitwise matrix multiplication 
		for i in range(len(result)):
			for j in range(26):
				# AND MULTIPLICATION
				mat_mult = check_block[j] and parity_matrix[j, i]
				# XOR ADDITION
				result[i] = (result[i] and not mat_mult) or (not result[i] and mat_mult)
		

		result = result.astype(int)
		result = (result).tolist()

		if (result == [1,1,1,1,0,1,1,0,0,0]):
			syndrome = "A"
			print("Syndrome is A at position: ", pos+26)
		
		if (result == [1,1,1,1,0,1,0,1,0,0]):
			syndrome = "B"
			print("Syndrome is B at position: ", pos+26)
		
		if (result == [1,0,0,1,0,1,1,1,0,0]):
			syndrome = "C"
			print("Syndrome is C at position: ", pos+26)
		
		if (result == [1,1,1,1,0,0,1,1,0,0]):
			syndrome = "C prime"
			print("Syndrome is C' at position: ", pos+26)
		
		if (result == [1,0,0,1,0,1,1,0,0,0]):
			syndrome = "D"
			print("Syndrome is D at position: ", pos+26)

		
		pos += 1
		if pos+26 > len(diff_stream):
			synced = True
