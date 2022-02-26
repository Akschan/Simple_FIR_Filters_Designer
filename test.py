import numpy as np
import matplotlib.pyplot as plt
dt = 0.001
t = np.arange(-1,1,dt)
w = 2 * np.pi * 250
f = np.sin(w*t)/(np.pi*t)
f = f * dt
#f = f + 2.8 * np.random.randn(len(t))
n = len(t)
fft = np.fft.fft(f,n)
PSD = np.log10(np.real(fft)*np.real(fft) + np.imag(fft)*np.imag(fft))
freqvect =(1/(dt*n)) * np.arange(n)
l = np.arange(1,np.floor(n/2),dtype='int')

fig, axs = plt.subplots(2)
axs[0].plot(t,f)
axs[1].plot(freqvect,PSD)
plt.show()
