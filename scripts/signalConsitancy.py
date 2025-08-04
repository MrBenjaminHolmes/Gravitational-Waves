file_name = 'PyCBC_T3_0.gwf'
channel_name = "H1:TEST-STRAIN"
start = 0
end = start + 128

get_file(file_name)

# Load and preprocess strain
ts = read_frame(file_name, channel_name, start, end)
ts = ts.highpass_fir(15, 512)
ts = resample_to_delta_t(ts, 1.0/2048).crop(2, 2)

# Compute PSD
psd = ts.psd(2)
psd = interpolate(psd, ts.delta_f)
psd = inverse_spectrum_truncation(psd, int(2 * ts.sample_rate), low_frequency_cutoff=15.0)

# Generate waveform
hp, _ = get_fd_waveform(approximant="IMRPhenomD",
                         mass1=32, mass2=32,
                         f_lower=20.0, delta_f=ts.delta_f)
hp.resize(len(ts.to_frequencyseries()))

# Matched filter
snr = matched_filter(hp, ts, psd=psd, low_frequency_cutoff=20)
snr = snr.crop(5, 4)

# Plot SNR
plt.plot(snr.sample_times, abs(snr.data))
plt.title("SNR")
plt.xlim(100, 120)
plt.show()

# Chi-squared test
nbins = 26
chisq = power_chisq(hp, ts, nbins, psd, low_frequency_cutoff=20.0)
chisq = chisq.crop(5, 4)
dof = nbins * 2 - 2
chisq /= dof 

plt.figure(figsize=[14, 4])
plt.plot(chisq.sample_times, chisq.data)
plt.grid()
plt.title("$\chi^2_r$ Consistency test")
plt.xlim(100, 120)
plt.ylabel('$\chi^2_r$')
plt.show()


nsnr_x = newsnr(abs(snr), chisq)

# Plot the new SNR timeseries
plt.figure(figsize=[14, 4])
plt.plot(snr.sample_times, nsnr_x, label='H1')
plt.title('NewSNR Timeseries')
plt.grid()
plt.xlim(100,120)
plt.ylim(0, 15)
plt.xlabel('Time (s)')
plt.ylabel('Re-weighted Signal-to-noise')
plt.show()
