import sys
import re
import random
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 4 or len(sys.argv) != int(sys.argv[2])+4 :
	print ("Program syntax: gen.py [duration in sec] [number of] [sinus | ripples | points] [freqence1:ampTab1 frequence2:ampTab2 ...]")
	exit()


frameRate = 44100
duration = float(sys.argv[1])
numOfWaves = int(sys.argv[2])
waveType = sys.argv[3]
x = np.arange(0, duration, 1/frameRate)
f = np.zeros(int(duration*frameRate))
for j in range(numOfWaves):
	wave = sys.argv[j+4]
	waveArgv = re.split(":", wave)
	freqBase = int(waveArgv[0])
	ampArray = ["1000"]
	if len(waveArgv) == 2: 
		ampArray = re.split(",", waveArgv[1])
	for i in range(len(ampArray)):
		freq = freqBase * (i+1)
		amp = int(ampArray[i])
		phase = random.randint(0, int(frameRate/freq))
		if waveType == "sinus":
			for i in range(int(duration*frameRate)):
				f[i] += amp*np.sin((i + phase) * freq  / frameRate * 2 * np.pi)
		if waveType == "ripples":
			a = 1
			for i in range(int(duration*frameRate)):
				f[i] += a*np.cosh((((i + phase) * freq / frameRate) % 1 - 0.5) * 2 * np.pi / a)
		if waveType == "points":
			for i in range(int(duration*frameRate)):
				if (i + phase) % int(frameRate/freq) == 0: f[i] += amp
	
	
plt.plot(x, f)	
plt.xlabel("time [s]")
plt.ylabel("amplitude")
plt.show()

outputFile = open("wave.in", "w")
print(frameRate, file=outputFile)
print(duration, file=outputFile)
print(*f, file=outputFile)
outputFile.close()
