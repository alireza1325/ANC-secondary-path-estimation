# ANC-secondary-path-estimation

The app is designed to tackle a key challenge in the development of active noise cancellation (ANC) systems, which is the presence of nonlinearities in the secondary path that can limit system performance. The secondary path can be modeled using a Finite Impulse Response (FIR) filter.

The app's main.cpp file implements an adaptive FIR filter to estimate the acoustic secondary path in an ANC system. It generates a known signal for 5 seconds through a loudspeaker and records it using a microphone. The FIR filter coefficients are updated in real-time using the portaudio library in C++.

To use the app, you need to first find the index of your input and output devices by compiling and running the pa_devs.c provided by portaudio. Then, you can put these indexes in the main.cpp file.

The app saves the input, recorded, and filtered sounds to WAV files, which can be used for further analysis or processing. 

## How to compile the code?

First you need to install the portaudio library in your computer using the information on this page: http://www.portaudio.com/docs/v19-doxydocs/compile_linux.html

Then run this command to compile the code:
```
g++ -o main main.cpp -lrt -lasound -ljack -lpthread -lportaudio -lsndfile
```
Make sure you have installed g++ compiler before in your device. You can also install the sndfile library for manafing the wav files using this command: 

```
apt-get install libsndfile-dev
```


