import wave
import struct 
import numpy as np
import cmath 
import math
import matplotlib.pyplot as plt
import pylab

contraOctave = {"C" : 32.70, "Cis" : 34.65, "D" : 36.71, "Dis" : 38.89,
"E" : 41.20, "F" : 43.65, "Fis" : 46.25, "G" : 49.00, "Gis" : 51.91,
"A" : 55.00, "Ais" : 58.27, "H" : 61.74} #tone : frequence

inf = 2000000000
fourierSection = 50 #ms
track = wave.open("Vivaldi-Sovente-il-sole.wav", mode="rb")
length = track.getnframes()
frameRate = track.getframerate()
fouriers = 300;

#reading file
samples = int(fouriers*frameRate*fourierSection/1000)
data = []
for sample in range(samples):
	waveData = track.readframes(1)
	probe = struct.unpack("<2h", waveData)
	data.append(probe[0]) #channel 0, quit channel 1

#fft
fourierStep = int(frameRate*fourierSection/1000)
toneTab = { i : [] for i in contraOctave.keys() }
for key, val in contraOctave.items():
	fourierDuration = 1/val
	fourierSize = int(frameRate*fourierDuration)
	for i in range(fouriers):
		freqRow = []
		dataRow = []
		for j in range(fourierSize):
			dataRow.append(data[i*fourierStep+j])
		dft = np.fft.rfft(dataRow)
		for j in range(len(dataRow)//2):
			freqRow.append(np.log2(abs(dft[j])+1))
		toneTab[key].append(freqRow)

freqTab = []
minimalToneFreqRowSize = int(frameRate*1/contraOctave["C"]//4)
print(minimalToneFreqRowSize)
for i in range(fouriers):
	freqRow = []
	for j in range(minimalToneFreqRowSize):
		for tone, toneFreqRowArray in toneTab.items():
			freqRow.append(toneFreqRowArray[i][j])
	freqTab.append(freqRow)
		
#spectrogram
specAxisY = []
for i in range(minimalToneFreqRowSize):
	for tone, freq in contraOctave.items():
		specAxisY.append(freq*i)
		
amp = np.array(freqTab).transpose()
plt.subplot(2, 1, 1)
plt.pcolormesh(np.arange(0, fouriers/20, 1/20), specAxisY, amp)
plt.xlim(5,15)
plt.ylim(32,11000)
plt.ylabel('frequency [Hz]')
plt.xlabel('time [s]')


#difference diagram
diffArray = []
for i in range(fouriers-1):
	sum = 0
	for j in range(len(specAxisY)):
		sum += abs(freqTab[i+1][j]-freqTab[i][j])/(freqTab[i+1][j]+freqTab[i][j]+1)
	diffArray.append(sum)
sortedDiffArray = sorted(diffArray)
topVal = sortedDiffArray[int(len(diffArray)*0.99)]
bottomVal = sortedDiffArray[int(len(diffArray)*0.01)]
plt.subplot(2, 1, 2)
plt.plot(np.arange(0, (fouriers-1)/20, 1/20), diffArray)
plt.xlim(5,15)
plt.ylim(bottomVal,topVal)
plt.xlabel('time [s]')

plt.show()
