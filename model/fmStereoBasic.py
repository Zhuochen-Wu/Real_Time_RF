import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math
from fmSupportLib import fmDemodArctan,estimatePSD
from fmPll import fmPll 


# the radio-frequency (RF) sampling rate
# this sampling rate is either configured on RF hardware
# or documented when a raw file with IQ samples is provided
rf_Fs = 2.4e6

# the cutoff frequency to extract the FM channel from raw IQ data
rf_Fc = 100e3

# the number of taps for the low-pass filter to extract the FM channel
# this default value for the width of the impulse response should be changed
# depending on some target objectives, like the width of the transition band
# and/or the minimum expected attenuation from the pass to the stop band
rf_taps = 151

# the decimation rate when reducing the front end sampling rate (i.e., RF)
# to a smaller samping rate at the intermediate frequency (IF) where
# the demodulated data will be split into the mono/stereo/radio data channels
rf_decim = 10

# audio sampling rate (we assume audio will be at 48 KSamples/sec)
audio_Fs = 240e3

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

block_size = 1024 * rf_decim * audio_decim * 2

def customlfilter(b, a, x):
    a = 1.0
    N = len(b)
    Yn = [0]*len(x)

    for n in range(len(x)):
    	for i in range(N):
    		# print(i)
    		if (n < i):
    			break
    		else:
    			Yn[n] += b[i]*x[n-i]

    return Yn


def bandpassFilter(Fb, Fe, Fs, num_taps):
    norm_center = ((Fe+Fb)/2)/(Fs/2)
    norm_pass = (Fe-Fb)/(Fs/2)

    h = [0] * (num_taps-1) 
    for i in range(num_taps-1):
        if (i == (num_taps - 1)/2):
            h[i] = norm_pass
        else:
            h[i] = norm_pass * ((math.sin(math.pi*(norm_pass/2)*(i-(num_taps-1)/2))) / (math.pi*(norm_pass/2)*(i-(num_taps-1)/2)))
    
        h[i] = h[i] * math.cos(i*math.pi*norm_center)
        h[i] = h[i] * pow(math.sin((i*math.pi)/num_taps), 2)

    return h


# custom lowpass filter
def lowPassFilter(Fc, Fs, Ntaps):
	norm_Fc = Fc/(Fs/2)
	h = [0]*(Ntaps-1)

	for i in range(Ntaps-1):
		if i == (Ntaps-1)/2:
			h[i] = norm_Fc
		else:
			h[i] = norm_Fc*(math.sin(math.pi*norm_Fc*(i-(Ntaps-1)/2))/(math.pi*norm_Fc*(i-(Ntaps-1)/2)))

		h[i] = h[i]*(math.sin(i*math.pi/Ntaps))**2

	return h


def plotSpectrum(x, Fs, type = 'FFT'):

	n = len(x)             # length of the signal
	df = Fs/n              # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)
	#XMag = abs(Xfreq)/n
		
	if type == 'DFT':
		Xfreq = DFT(x)
    

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = abs(Xfreq)/n
	XMag = XMag[0:int(n/2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs/2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
		title='Frequency domain plot')
	# fig.savefig("freq.png")
	plt.show()

def plotTime(x, time):

    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()


	

if __name__ == "__main__":
     in_fname = "../data/test1.raw"
     iq_data = np.fromfile(in_fname, dtype='float32')
     print("Read raw RF data from \"" + in_fname + "\" in float32 format")
     #Low Pass filter Coeff's to extract fm Channel 
     
     #rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
     rf_coeff = lowPassFilter(rf_Fc, rf_Fs, rf_taps)    
     i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0:block_size:2])
     q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1:block_size:2])
     
     i_ds = i_filt[::rf_decim]
     q_ds = q_filt[::rf_decim]
     
     fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
     
    #  nyquist = audio_Fs//2
    #  norm_cutoff = audio_fc/nyquist

    #for pilot tone it should be 18.5e3 and 19.5e3
     

     #bandPassCoeff_pilot = signal.firwin(rf_taps, [18.5e3/(audio_Fs/2), 19.5e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
     #print(bandPassCoeff_pilot)
     print("\n our coeff: \n")
     # bpf_recovery = signal.lfilter(bandPassCoeff_pilot, 1.0, fm_demod)
     #####
     ourBandPassCoeff = bandpassFilter(18.5e3, 19.5e3, audio_Fs, rf_taps)
     bpf_recovery = signal.lfilter(ourBandPassCoeff, 1.0, fm_demod)     

     # print(ourBandPassCoeff)

     recovery_pll = fmPll(bpf_recovery, 19e3, 240e3,2)

     print("Completed Recovery BPF Convolution")

     #bandPassCoeff_stereo = signal.firwin(rf_taps, [22e3/(audio_Fs/2), 54e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
     #bpf_recovery_stereo = signal.lfilter(bandPassCoeff_stereo, 1.0, fm_demod)
     ourBandPassCoeff_stereo = bandpassFilter(22e3, 54e3, audio_Fs, rf_taps)
     bpf_recovery_stereo = signal.lfilter(ourBandPassCoeff_stereo, 1.0, fm_demod)
     print("Completed Stereo Channel Convolution Block 1")
     
     fig, (ax0,ax1,ax2) = plt.subplots(nrows=3)
     fig.subplots_adjust(hspace = 1.0)

     ax0.psd(bpf_recovery, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
     ax0.set_ylabel('PSD (db/Hz)')
     ax0.set_title('Extracted Stereo Pilot Tone')

     ax1.psd(bpf_recovery_stereo, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
     ax1.set_ylabel('PSD (db/Hz)')
     ax1.set_title('Extracted Stereo Channel')

     
     #MIXINGG
     #fmpll(x,pllfreq,NCO Freq,FS)
     
     #plotSpectrum(recovery_pll,240e3)

     mixed = np.multiply(recovery_pll[0:len(bpf_recovery_stereo):1],bpf_recovery_stereo)
    
    
    #  for x in range(len(recovery_pll)):
    #      mixed recovery_pll[x]*bpf_recovery_stereo[x]
        
        
     # stereo_coeff = signal.firwin(rf_taps, 16e3/(audio_Fs/2), window=('hann'))
     
     ##CHANGE TO AUDIO COEFFS NOT RF COEFFICENTS 
     
     stereo_coeff = lowPassFilter(audio_Fc, audio_Fs, audio_taps)
     stereo_filt = signal.lfilter(stereo_coeff, 1.0, mixed)  
     #Decimation
     stereo_data = stereo_filt[::5] 
    
     # PSD after Downconversion and filtering the Stereo Channel
     ax2.psd(stereo_filt, NFFT=512, Fs=audio_Fs/1e3)
     ax2.set_ylabel('PSD (db/Hz)')
     ax2.set_title('Downconversion and PLL mixed')

     #lotSpectrum(recovery_pll[0:len(recovery_pll-1):1],240e3)

    


    # mono audio
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	fig.subplots_adjust(hspace = 1.0)

	# PSD after FM demodulation
	ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax0.set_ylabel('PSD (db/Hz)')
	ax0.set_title('Demodulated FM')

	# coefficients for the filter to extract mono audio
	# audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann')) # to be updated by you during in-lab
	audio_coeff = FIR(audio_Fc, audio_Fs, audio_taps)
	# extract the mono audtio data through filtering

	#Ask if we should just due the decimation here so we do fm_demod[::5]
	#audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod) # to be updated by you during in-lab
	
	
	# audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod[::audio_decim]) # to be updated by you during in-lab
	audio_filt =  convolution(audio_coeff, fm_demod[::audio_decim])
	# you should uncomment the plots below once you have processed the data

	#PSD after extracting mono audio
	ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
	ax1.set_ylabel('PSD (db/Hz)')
	ax1.set_title('Extracted Mono')

	# downsample audio data
	#audio_data = audio_filt[::audio_decim] # to be updated by you during in-lab
	audio_data = audio_filt

	#PSD after decimating mono audio
	ax2.psd(audio_data, NFFT=512, Fs=audio_Fs/1e3)
	ax2.set_ylabel('PSD (db/Hz)')
	ax2.set_title('Mono Audio')

	# save PSD plots
	fig.savefig("../data/fmMonoBasic.png")

    # combiner
    ### left and right
    left_channel = np.array(len(audio_data))
    right_channel = np.array(len(audio_data))
    for i in range(len(audio_data)):
        left_channel[i] = (stereo_data[i] + audio_data[i])/2
        right_channel[i] = (audio_data[i] - stereo_data[i])/2
    


    plt.show()