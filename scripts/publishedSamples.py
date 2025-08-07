label = 'GW150914'
! wget --no-verbose https://dcc.ligo.org/LIGO-P1800370/public/{label}_GWTC-1.hdf5
posterior_file = './'+label+'_GWTC-1.hdf5'
posterior = h5py.File(posterior_file, 'r')
print('This file contains four datasets: ',posterior.keys())
print("Properties: ",posterior['Overall_posterior'].dtype.names)

samples=pd.DataFrame.from_records(numpy.array(posterior['Overall_posterior']))
samples
corner.corner(samples.values,labels=['costhetajn',
                                'distance [Mpc]',
                                'ra',
                                'dec',
                                'mass1 [Msun]',
                                'mass2 [Msun]',
                                'spin1',
                                'spin2',
                                'costilt1',
                                'costilt2']);
