def get_file(fname):
    url = f"https://raw.githubusercontent.com/MrBenjaminHolmes/Gravitational-Waves/main/SampleData/2.22/{fname}"
    urllib.request.urlretrieve(url, fname)
    print(f'Getting: {url}')

files = ['PyCBC_T2_0.gwf', 'PyCBC_T2_1.gwf', 'PyCBC_T2_2.gwf']
channel_name = "H1:TEST-STRAIN"
start = 0
end = start + 128
conditionedData = []
fig, axs = plt.subplots(len(files), 1, figsize=(10, 7), sharex=True)

for i, file_name in enumerate(files):
    get_file(file_name)
    ts = read_frame(file_name, channel_name, start, end)
    ts  = highpass(ts ,15)
    ts  = resample_to_delta_t(ts ,1/2048)
    ts = ts.crop(2,2) 
    conditionedData.append([file_name,ts])
    axs[i].plot(ts.sample_times, ts.data)
    axs[i].set_title(f"{file_name}")
    axs[i].grid(False)
    axs[i].set_ylabel("Strain")

axs[-1].set_xlabel("Time (s)")
plt.tight_layout()

plt.show()

masses = [15, 30, 45]

for i, strain_data in enumerate(conditionedData):
    fig, axs = plt.subplots(1, len(masses), figsize=(15, 5), sharey=True)
    fig.suptitle(strain_data[0])
    
    print(f"------------- {strain_data[0]} -------------")
    
    for j, mass in enumerate(masses):
        # Create PSD
        psd = strain_data[1].psd(4)
        psd = interpolate(psd,  strain_data[1].delta_f)
        psd = inverse_spectrum_truncation(psd, int(4 *  strain_data[1].sample_rate),
                                          low_frequency_cutoff=15)

        # Create Event Template
        hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                                 mass1=mass,
                                 mass2=mass,
                                 delta_t= strain_data[1].delta_t,
                                 f_lower=20)

        # Resize waveform vector to match data length
        hp.resize(len( strain_data[1]))

        # Matched filter
        snr = matched_filter(hp,strain_data[1],
                             psd=psd, low_frequency_cutoff=20)
        snr = snr.crop(4 + 4, 4)

        peak = abs(snr).numpy().argmax()
        snrp = snr[peak]
        time = snr.sample_times[peak]

        axs[j].plot(snr.sample_times, abs(snr))
        axs[j].set_title(f'Mass = {mass}')
        axs[j].set_xlabel('Time (s)')
        if j == 0:
            axs[j].set_ylabel('SNR')

        print(f"Signal at {time}s with SNR {abs(snrp)} using mass {mass}")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
