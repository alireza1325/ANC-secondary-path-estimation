/** @file main2.cpp
    @ingroup examples_src
    @brief just playback and record what it heard
    @author Alireza
*/
/*
 * $Id$
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 */

// This app, play a sound and record it back using microphone to estimate the secondary path in ANC system using FIR filter.

#include <stdio.h>
#include <math.h>
#include "portaudio.h"
#include "sndfile.h"
#include <vector>
#include <iostream>
#include <numeric> // for std::inner_product


#define SAMPLE_RATE         (44100)
#define PA_SAMPLE_TYPE      paFloat32
#define FRAMES_PER_BUFFER   (256)
#define NUM_CHANNELS    (1)
#define NUM_SECONDS     (5)
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE   (1)

typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"

float systemDelay = 0.041;  //0.041 s
int sampleDelay = systemDelay*SAMPLE_RATE;

typedef struct
{
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
    SAMPLE      *playSamples;
    PaTime streamTime;
    PaTime inputLatency;
    PaTime outputLatency;

}
paTestData;


void toWav(SAMPLE* data,const char* filename)
{
    //const char* filename = "output.wav";
    const int sample_rate = SAMPLE_RATE;
    const float duration = NUM_SECONDS;
    const int num_samples = sample_rate * duration;

    SF_INFO sfinfo;
    sfinfo.channels = NUM_CHANNELS;
    sfinfo.samplerate = sample_rate;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;

    SNDFILE* outfile = sf_open(filename, SFM_WRITE, &sfinfo);
    if (!outfile) {
        printf("Error opening WAV file for writing: %s\n", sf_strerror(NULL));
    }

    sf_count_t framesWritten = sf_write_float(outfile, data, num_samples);
    if (framesWritten != num_samples) {
        printf("Error writing WAV file: %s\n", sf_strerror(outfile));
    }
    sf_close(outfile);

}

// FxLMS function for secondary path estimation

// define constants
#define NUM_TAPS (256)
#define MU (0.008)

// Adaptive FIR filter implementation
class AdaptiveFIRFilter {
public:
  AdaptiveFIRFilter(int filter_length, double step_size)
      : filter_length_(filter_length), step_size_(step_size) {
    filter_.resize(filter_length_, 0.0);
    input_history_.resize(filter_length_, 0.0);
  }

  void ApplyLMS(float d,float input) {

    // update the noise history by the new input 

    //input_history_.insert(input_history_.begin(), input);
    //input_history_.pop_back();
    
    for (int i = filter_length_ - 1; i > 0; i--)
    {
        input_history_[i] = input_history_[i-1];
    }
    input_history_[0] = input;

    float output = 0.0;
    for (int i = 0; i < filter_length_; ++i) {
      output += filter_[i] * input_history_[i];
    }
    float error = d - output;
    for (int i = 0; i < filter_length_; ++i) {
      filter_[i] += step_size_ * error * input_history_[i];
    }
  }

  std::vector<float> applyFIRFilter(const std::vector<float>& signal)
  {
    // Initialize output vector with same length as input signal
    std::vector<float> filteredSignal(signal.size(), 0.0);

    // Apply filter
    for (int n = 0; n < signal.size(); n++) {
        // Calculate output at current time step
        float y_n = 0;
        for (int k = 0; k < filter_length_; k++) {
            if (n - k >= 0) {
                y_n += filter_[k] * signal[n - k];
            }
        }
        filteredSignal[n] = y_n;
    }
    return filteredSignal;
  }

  std::vector<float> GetFilter(){
    return filter_;
  }

private:
  int filter_length_;
  float step_size_;
  std::vector<float> filter_;
  std::vector<float> input_history_;
};


AdaptiveFIRFilter filt(NUM_TAPS,MU);



static int gNumNoInputs = 0;
/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int SWCallback( const void *inputBuffer, void *outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void *userData )
{

    paTestData *data = (paTestData*)userData;
    //const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr2 = &data->playSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;

    SAMPLE *out = (SAMPLE*)outputBuffer;
    const SAMPLE *in = (const SAMPLE*)inputBuffer;
    
    // get timing information
    data->streamTime = timeInfo->currentTime; // Get current time
    data->inputLatency = timeInfo->inputBufferAdcTime - timeInfo->currentTime; // Compute input latency
    data->outputLatency = timeInfo->outputBufferDacTime - timeInfo->currentTime; // Compute output latency

    // end timing information

    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

    (void) timeInfo; /* Prevent unused variable warnings. */
    (void) statusFlags;
    (void) userData;

    

    if( framesLeft < framesPerBuffer )
    {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else
    {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }

    if( inputBuffer == NULL )
    {
        for( i=0; i<framesToCalc; i++ )
        {
            *out++ = 0;  /* left - silent */
            *wptr++ = SAMPLE_SILENCE;  /* left */
            *wptr2++ = SAMPLE_SILENCE;

            if( NUM_CHANNELS == 2 )
            {
                *out++ = 0;  /* right - silent */
                *wptr++ = SAMPLE_SILENCE;
                *wptr2++ = SAMPLE_SILENCE;
            }
             
        }
        gNumNoInputs += 1;
    }
    else
    {
        for( i=0; i<framesToCalc; i++ )
        {
            SAMPLE sample = *in++; /* MONO input */
            *out++ = *wptr2++;
            *wptr++ = sample;  /* left */

            if(data->frameIndex>sampleDelay) filt.ApplyLMS(sample,data->playSamples[(data->frameIndex-sampleDelay)* NUM_CHANNELS]);

            if( NUM_CHANNELS == 2 )
            {
                SAMPLE sample = *in++;
                *out++ = *wptr2++;       /* right */
                *wptr++ = sample;  /* right */
                
                if(data->frameIndex>sampleDelay) filt.ApplyLMS(sample,data->playSamples[(data->frameIndex-sampleDelay)* NUM_CHANNELS]);
            } 
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

/*******************************************************************/
int main(void);
int main(void)
{


    
    // read the file you want to play
    SF_INFO sfinfo;
    sfinfo.format = 0;
    SNDFILE* file = sf_open("10m_white_noise_44100Hz.wav", SFM_READ, &sfinfo);
    if (!file) {
        std::cout << "Error: failed to open audio file" << std::endl;
        return 1;
    }

    // Read the audio data
    std::vector<_Float32> audioData(sfinfo.frames);
    sf_read_float(file, audioData.data(), sfinfo.frames);

    // Close the file
    sf_close(file);

    std::cout << "Total numer of frames: "<<audioData.size()<< ", Which is "<<  audioData.size()/SAMPLE_RATE <<" s" <<std::endl;


    // initialized the portaudio library
    PaStreamParameters inputParameters, outputParameters;
    PaStream *stream;
    PaError err;
    paTestData          data;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;

    
    
    data.maxFrameIndex = totalFrames = (float)(NUM_SECONDS +systemDelay)* SAMPLE_RATE; /* Record for a few seconds. */
    data.frameIndex = 0;
    numSamples = totalFrames * NUM_CHANNELS;
    numBytes = numSamples * sizeof(SAMPLE);
    data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
    data.playSamples = (SAMPLE *) malloc( numBytes );
    

    std::vector<float> sw(NUM_TAPS,0.0);
    std::vector<float> filtredSignal(numSamples);
    std::vector<float> inputSignal(numSamples);
    SAMPLE filtredSignals[filtredSignal.size()];

    
    if( data.recordedSamples == NULL )
    {
        printf("Could not allocate record array.\n");
        goto error;
    }
    for( int i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;
    for( int i=0; i<numSamples; i++ )
    {
        if (i>systemDelay*SAMPLE_RATE) data.playSamples[i] = audioData[i];
        else data.playSamples[i] = 0;
    } 
    

    err = Pa_Initialize();
    if( err != paNoError ) goto error;

    inputParameters.device = 7; /* set input device */ 
    if (inputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default input device.\n");
        goto error;
    }
    inputParameters.channelCount = NUM_CHANNELS;       /* mono input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    outputParameters.device = 6; /* set output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto error;
    }
    outputParameters.channelCount = NUM_CHANNELS;       /* mono output */
    outputParameters.sampleFormat = PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;

    err = Pa_OpenStream(
              &stream,
              &inputParameters,
              &outputParameters,
              SAMPLE_RATE,
              FRAMES_PER_BUFFER,
              0, /* paClipOff, */  /* we won't output out of range samples so don't bother clipping them */
              SWCallback,
              &data);
    if( err != paNoError ) goto error;

    err = Pa_StartStream( stream );
    if( err != paNoError ) goto error;


    while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
    {
        Pa_Sleep(1000);
        printf("index = %d, streamTime = %f, inputLatency = %f, outputLatency = %f\n", data.frameIndex,data.streamTime,data.inputLatency,data.outputLatency); 
        fflush(stdout);

    }
    if( err < 0 ) goto error;
    

    //printf("Hit ENTER to stop program.\n");
    //getchar();
    err = Pa_CloseStream( stream );
    if( err != paNoError ) goto error;

    

#if WRITE_TO_FILE
    {
        toWav(data.recordedSamples,"output.wav");
        toWav(data.playSamples,"input.wav");
        printf("Wrote data to 'output.awv'\n");
    }
#endif

    for(int i=0;i<numSamples;i++)
    {
        inputSignal[i] = (float)data.playSamples[i];
    }

    sw = filt.GetFilter();
    printf("Secondary path weights\n");
    for(int i =0; i<NUM_TAPS;i++) printf("%f\t",sw[i]);
    printf("\n");
    
    filtredSignal = filt.applyFIRFilter(inputSignal);
    
    for(int i=0;i<numSamples;i++) filtredSignals[i] = filtredSignal[i];

    toWav(filtredSignals,"filteredSignal.wav");

    printf("Finished. gNumNoInputs = %d\n", gNumNoInputs );
    Pa_Terminate();


    if( data.recordedSamples )       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( data.playSamples )       /* Sure it is NULL or valid. */
        free( data.playSamples );
    return 0;

error:
    Pa_Terminate();
    if( data.recordedSamples)       /* Sure it is NULL or valid. */
        free( data.recordedSamples );
    if( data.playSamples )       /* Sure it is NULL or valid. */
        free( data.playSamples );
        
    fprintf( stderr, "An error occurred while using the portaudio stream\n" );
    fprintf( stderr, "Error number: %d\n", err );
    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
    return -1;
}
