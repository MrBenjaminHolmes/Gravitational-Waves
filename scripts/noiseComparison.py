#--GW190412 Data-#
gw19 = event_gps('GW190412')
gw19Data = TimeSeries.fetch_open_data('L1', int(gw19)-512, int(gw19)+512, cache=True)
gw19ASD = gw19Data.asd(fftlength=4, method="median")

#--GW150914 Data-#
gw15 = event_gps('GW150914')
gw15Data = TimeSeries.fetch_open_data('L1', int(gw15)-512, int(gw15)+512, cache=True)
gw15ASD = gw15Data.asd(fftlength=4, method="median")

plt.loglog(gw19ASD, label = 'GW190412')
plt.loglog(gw15ASD,label = 'GW150914')
plt.legend()

plt.show()
