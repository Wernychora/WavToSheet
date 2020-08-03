import wave
import struct 
import numpy as np
import cmath 
import math
import matplotlib.pyplot as plt
import pylab

fourierSection = 50 #ms
probeSection = 200 #ms
track = wave.open("Vivaldi-Sovente-il-sole.wav", mode="rb")
length = track.getnframes()
frameRate = track.getframerate()
fourierFrames = int(frameRate*fourierSection/1000)
fouriers = 1000;
freqTab = []
for i in range(fouriers):
	data = []
	freqRow = []
	for j in range(fourierFrames):
		waveData = track.readframes(1)
		probe = struct.unpack("<2h", waveData)
		#print(probe[0], probe[1])
		data.append(probe[0])
	#print("----------------------------")
	#print("data size", len(data))
	freq = np.fft.fftfreq(len(data), d=1./frameRate)
	dft = np.fft.rfft(data)
	for j in range(len(data)//2):
		#print(freq[j], ":", cmath.polar(dft[j])[0])
		freqRow.append(cmath.polar(dft[j])[0])
	#print("/////////////////////////////")
	freqTab.append(freqRow)

amp = np.array(freqTab).transpose()
plt.pcolormesh(np.arange(0, fouriers/20, 1/20), np.arange(0,22040,20), amp)
plt.ylabel('frequency [Hz]')
plt.xlabel('time [s]')
plt.show()



    
