import sys
import random
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
	print ("Program syntax: gen.py [duration in sec] [number of] [sinus | ripples | points]")
	exit()

frameRate = 44100
duration = float(sys.argv[1])
n = int(sys.argv[2])
waveType = sys.argv[3]
x = np.arange(0, duration, 1/frameRate)
f = np.zeros(int(duration*frameRate))
for j in range(n):
	freq = random.randint(64,1024)
	amp = random.randint(300, 1000)
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

print(frameRate)
print(duration)
print(*f)
