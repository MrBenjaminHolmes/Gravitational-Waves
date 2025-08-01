gps = event_gps("GW170817")
ldata = TimeSeries.fetch_open_data("L1", gps - 300, gps + 300, cache=True)
specgram = ldata.spectrogram2(fftlength=4, overlap=2, window='hann') ** (0.5)
plot = specgram.plot();

ax = plot.gca()
ax.set_yscale('log')
ax.set_ylim(10, 1400)
ax.colorbar(
    clim=(1e-24, 1e-20),
    norm="log",
    label=r"Strain noise [$1/\sqrt{\mathrm{Hz}}$]",
)
