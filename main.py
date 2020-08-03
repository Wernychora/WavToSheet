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
		freqRow.append(np.log2(cmath.polar(dft[j])[0]+1))
	#print("/////////////////////////////")
	freqTab.append(freqRow)
	
#spectrogram
amp = np.array(freqTab).transpose()
plt.subplot(2, 1, 1)
plt.pcolormesh(np.arange(0, fouriers/20, 1/20), np.arange(0,22040,20), amp)
plt.ylabel('frequency [Hz]')
plt.xlabel('time [s]')


#difference diagram
diffTable = []
for i in range(fouriers-1):
	sum = 0
	for j in range(len(data)//2):
		sum += abs(freqTab[i+1][j]-freqTab[i][j])/(freqTab[i+1][j]+freqTab[i][j]+1)
	diffTable.append(sum)
plt.subplot(2, 1, 2)
plt.plot(np.arange(0, (fouriers-1)/20, 1/20), diffTable)

plt.show()
	


    
