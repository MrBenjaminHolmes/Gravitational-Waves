#Get Merger Raw Signal
merger = Merger("GW150914")
strain = merger.strain("H1")

#Remove Low Frequencies with a Highpass
strain = highpass(strain,15)
#Resample data
strain = resample_to_delta_t(strain,1/2048)
#Crop Start and End Peaks
strain = strain.crop(2, 2)

plt.plot(strain.sample_times, strain)
plt.title("Raw Data")
plt.grid(False)
plt.xlabel('Time (s)')
plt.show()

#Create PSD
psd = strain.psd(4)
psd = interpolate(psd, strain.delta_f)
psd = inverse_spectrum_truncation(psd, int(4 * strain.sample_rate),
                                  low_frequency_cutoff=15)

#Create Event Template
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                     mass1=36,
                     mass2=36,
                     delta_t=strain.delta_t,
                     f_lower=20)

# Resize the Waveform Vector to match our data
hp.resize(len(strain))

snr = matched_filter(hp, strain,
                     psd=psd, low_frequency_cutoff=20)

snr = snr.crop(4 + 4, 4)
plt.plot(snr.sample_times,abs(snr))
plt.title("SNR")
plt.grid(False)
plt.xlabel("Time (s)")
plt.ylabel("Signal-to-Noise")
plt.show()

peak = abs(snr).numpy().argmax()
snrp = snr[peak]
time = snr.sample_times[peak]

print("We found a signal at {}s with SNR {}".format(time, abs(snrp)))

dt = time-strain.start_time
aligned = hp.cyclic_time_shift(dt)

#Scale so SNR=1 in the Data
aligned /= sigma(aligned, psd=psd, low_frequency_cutoff=20.0)

aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
aligned.start_time = strain.start_time


white_data = (strain.to_frequencyseries() / psd**0.5).to_timeseries()
white_template = (aligned.to_frequencyseries() / psd**0.5).to_timeseries()

#Highpass and Lowpass on Both
white_data = white_data.highpass_fir(30, 512).lowpass_fir(300, 512)
white_template = white_template.highpass_fir(30, 512).lowpass_fir(300, 512)

#Select the time around the merger
white_data = white_data.time_slice(merger.time-.2, merger.time+.1)
white_template = white_template.time_slice(merger.time-.2, merger.time+.1)


plt.figure(figsize=[15, 3])
plt.plot(white_data.sample_times, white_data, label="Data")
plt.plot(white_template.sample_times, white_template, label="Template")
plt.grid(False)
plt.title("Template VS Data")
plt.legend()
plt.show()

subtracted = strain-aligned


for data, title in [(strain, 'Original H1 Data'),
                    (subtracted, 'Signal Subtracted from H1 Data')]:
                      time,frequency,power = data.whiten(4, 4).qtransform(.001, logfsteps=100, qrange=(8, 8), frange=(20, 512))
                      plt.figure(figsize=[15, 3])
                      plt.pcolormesh(time,frequency,power**0.5, vmin=1, vmax=6, shading='auto')
                      plt.title(title)
                      plt.yscale('log')
                      plt.xlabel('Time (s)')
                      plt.ylabel('Frequency (Hz)')
                      plt.xlim(merger.time - 2, merger.time + 1)
                      plt.show()
