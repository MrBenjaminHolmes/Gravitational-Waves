gps = event_gps("GW170814")
segment = (int(gps) - 5, int(gps) + 2)
hData = TimeSeries.fetch_open_data('H1', *segment, verbose=True, cache=True)
lData = TimeSeries.fetch_open_data('L1', *segment, verbose=True, cache=True)
vData = TimeSeries.fetch_open_data('V1', *segment, verbose=True, cache=True)

hq = hData.q_transform(frange=(30, 400), qrange=(4, 20), outseg=(gps-0.2, gps+0.2))
lq = lData.q_transform(frange=(30, 400), qrange=(4, 20), outseg=(gps-0.2, gps+0.2))
vq = vData.q_transform(frange=(30, 400), qrange=(4, 20), outseg=(gps-0.2, gps+0.2))

hq_plot = hq.abs()**0.5
lq_plot = lq.abs()**0.5
vq_plot = vq.abs()**0.5

#-----------Graphing Data-----------#
fig, ax = plt.subplots(ncols=3, figsize=(15, 5), sharey=True)

xt_h = hq_plot.times.value - gps
xt_l = lq_plot.times.value - gps
xt_v = vq_plot.times.value - gps

im0 = ax[0].pcolormesh(xt_h, hq_plot.frequencies.value, hq_plot.value.T, vmin=1, vmax=5)
im1 = ax[1].pcolormesh(xt_l, lq_plot.frequencies.value, lq_plot.value.T, vmin=1, vmax=5)
im2 = ax[2].pcolormesh(xt_v, vq_plot.frequencies.value, vq_plot.value.T, vmin=1, vmax=5)

# Axis formatting
ax[0].set_title('H1')
ax[1].set_title('L1')
ax[2].set_title('V1')

ax[1].set_xlabel("Time [s]")

for a in ax:
    a.set_yscale('log')

ax[0].set_ylabel('Frequency [Hz]')


cbar = fig.colorbar(im2, ax=ax[2], pad=0.02)
cbar.set_label('Normalized energy')

fig.subplots_adjust(wspace=0.05)
plt.suptitle("Q-Transforms of GW170814 (H1, L1, V1)", fontsize=14)
plt.show()
#---------------------------------------#
