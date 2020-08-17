import wave
import struct 
import numpy as np
import cmath 
import math
import matplotlib.pyplot as plt
import pylab
import sys
from colorsys import hls_to_rgb

def read_file(data):
	samples = int(fouriers*frameRate*fourierSection/1000)
	for sample in range(samples):
		waveData = track.readframes(1)
		probe = struct.unpack("<2h", waveData)
		data.append(probe[0]) #channel 0, quit channel 1

def fft(freqTab, dft2Tab, toneTab, argData, argToneData):
	fourierStep = int(frameRate*fourierSection/1000)
	for key, val in contraOctave.items():
		fourierDuration = 1/val
		fourierSize = int(frameRate*fourierDuration)
		for i in range(fouriers):
			freqRow = []
			dataRow = []
			argRow = []
			window = np.hanning(fourierSize)
			for j in range(fourierSize):
				dataRow.append(data[i*fourierStep+j])
			dataRow = dataRow * window
			dft = np.fft.rfft(dataRow)
			for j in range(len(dataRow)//2):
				freqRow.append(np.log2(abs(dft[j])+1))
				argRow.append(np.angle(dft[j]))
			toneTab[key].append(freqRow)
			argToneData[key].append(argRow)
			
	minimalToneFreqRowSize = int(frameRate*1/contraOctave["C"]//4)

	for i in range(fouriers):
		freqRow = []
		argRow = []
		for j in range(minimalToneFreqRowSize):
			for tone, toneFreqRowArray in toneTab.items():
				freqRow.append(toneFreqRowArray[i][j])
			for tone, argToneArray in argToneData.items():
				argRow.append(argToneArray[i][j])
		freqTab.append(freqRow)
		argData.append(argRow)
	
	#fft2
	for i in range(fouriers):
		fftCol = []
		dft2 = np.fft.rfft(freqTab[i])
		dft2Tab.append(dft2)
	
	return minimalToneFreqRowSize

def diffPlot(array, flagArray, title, nr, bottomMargin, topMargin, differentiator = "null"):
	sortedArray = sorted(array)
	topVal = sortedArray[int(len(array)*topMargin)]
	bottomVal = sortedArray[int(len(array)*bottomMargin)]
	
	plt.subplot(2, 2, nr)
	plt.plot(np.arange(0, (fouriers-1)/20, 1/20), array)
	if differentiator == "zeroPoint":
		for x in range(len(array)):
			if array[x] < topVal/20 and array[x] > bottomVal/20:
				plt.plot(x/20, 0, marker='o', markersize=3, color="red")
				flagArray[x] = True
	if differentiator == "negativePoint":
		for x in range(len(array)):
			if array[x] < bottomVal/3:
				plt.plot(x/20, 0, marker='o', markersize=3, color="red")
				flagArray[x] = True
	plt.title(title)
	plt.xlim(0,15)
	plt.ylim(bottomVal,topVal)
	plt.xlabel("time [s]")

def colorize(r, arg):
    h = ((arg + np.pi) / (2 * np.pi) + 0.5) % 1.0
    l = r/np.amax(r)
    s = 0.8
    colorMatrix = np.vectorize(hls_to_rgb)(h,l,s)
    colorMatrix = np.array(colorMatrix)
    colorMatrix = colorMatrix.swapaxes(0,2) 
    return colorMatrix
    
def display(freqTab, dft2Tab, toneFourierSamples, argData, mode):
	#spectrogram data
	specAxisY = []
	for i in range(toneFourierSamples):
		for tone, freq in contraOctave.items():
			specAxisY.append(freq*i)		
	amp = np.array(freqTab)
	arg = np.array(argData)
	
	if mode == "fftComplex":
		plt.imshow(colorize(amp, arg), interpolation="none", extent=[0,15,0,11000], aspect="auto")
		plt.show()
		return

	#diffPlots
	diffArray = []
	diffFlags = []
	if mode == "diffPlots":
		#difference plot
		for i in range(fouriers-1):
			sum = 0
			for j in range(int(len(specAxisY)/10), len(specAxisY)):
				sum += abs(freqTab[i+1][j]-freqTab[i][j])/(freqTab[i+1][j]+freqTab[i][j]+1)
			diffArray.append(sum)	
		diffPlot(diffArray, diffFlags, "Difference diagram", 2, 0.01, 0.99)
		
		#diff derivative plot
		diffPrim = np.gradient(diffArray, edge_order = 2)
		diffPrimFlags = [False] * len(diffArray)
		diffPlot(diffPrim, diffPrimFlags, "Derivative diagram", 3, 0.01, 0.99, "zeroPoint")
		
		#diff second derivative plot
		diffBis = np.gradient(diffPrim, edge_order = 2)
		diffBisFlags = [False] * len(diffArray)
		diffPlot(diffBis, diffBisFlags, "Second derivative diagram", 4, 0.02, 0.98, "negativePoint")
		
	#spectrogram display
	if mode == "diffPlots": plt.subplot(2, 2, 1)
	if mode == "spectrograms": plt.subplot(2, 1, 1)
	plt.pcolormesh(np.arange(0, fouriers/20, 1/20), specAxisY, amp, shading="nearest")
	if mode == "diffPlots":
		for x in range(2, len(diffArray)-2):
			if diffPrimFlags[x] and (diffBisFlags[x-2] or diffBisFlags[x-1] or diffBisFlags[x] or diffBisFlags[x+1]):
				plt.axvline(x/20, 0, 10000, label='pyplot vertical line', color="red")
	plt.title("Spectrogram")
	plt.xlim(0,15)
	plt.ylim(0,11000)
	plt.ylabel("frequency [Hz]")
	plt.xlabel("time [s]")
	
	#fft2 spectrogram
	if mode == "spectrograms":
		dft2Data = np.array(dft2Tab).transpose()
		plt.subplot(2, 1, 2)
		plt.pcolormesh(np.arange(0, fouriers/20, 1/20), np.arange(0, len(specAxisY)), amp, shading="nearest")
		plt.title("DFT2 spectrogram")
		plt.xlim(0,15)
		plt.ylim(0,len(specAxisY))
		#plt.ylabel("frequency [Hz]")
		plt.xlabel("time [s]")
		
	plt.show()
	
#####
	
if len(sys.argv) != 3:
	print ("Program syntax: main.py [name of *.wav file] [spectrograms | diffPlots | fftComplex]")
	exit()
	
mode = sys.argv[2]

contraOctave = {"C" : 32.70, "Cis" : 34.65, "D" : 36.71, "Dis" : 38.89,
"E" : 41.20, "F" : 43.65, "Fis" : 46.25, "G" : 49.00, "Gis" : 51.91,
"A" : 55.00, "Ais" : 58.27, "H" : 61.74} #tone : frequence

inf = 2000000000
fourierSection = 50 #ms
track = wave.open(sys.argv[1], mode="rb")
length = track.getnframes()
frameRate = track.getframerate()
fouriers = 300;

#reading file
data = []
read_file(data)

#fft
freqTab = []
dft2Tab = []
argData = []
argToneData = { i : [] for i in contraOctave.keys() }
toneTab = { i : [] for i in contraOctave.keys() }
toneFourierSamples = fft(freqTab, dft2Tab, toneTab, argData, argToneData)

#spectrogram + difference plot (modes: "spectrograms", "diffPlots", "fftComplex")
display(freqTab, dft2Tab, toneFourierSamples, argData, mode)
