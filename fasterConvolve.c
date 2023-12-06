#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

typedef struct {
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short audio_format;
    short num_channels;
    int sample_rate;
    int byte_rate;
    short block_align;
    short bits_per_sample;
} WavHeader;

//Code inspired from the file in Course Documents in the Audio FIle Formats Test Tone Sample Code C file
void fwriteIntLSB(int value, FILE *file) {
    unsigned char buffer[4];
    buffer[0] = (value & 0x000000FF);
    buffer[1] = (value & 0x0000FF00) >> 8;
    buffer[2] = (value & 0x00FF0000) >> 16;
    buffer[3] = (value & 0xFF000000) >> 24;
    fwrite(buffer, 4, 1, file);
}
//Code inspired from the file in Course Documents in the Audio FIle Formats Test Tone Sample Code C file
void fwriteShortLSB(short value, FILE *file) {
    unsigned char buffer[2];
    buffer[0] = (value & 0x00FF);
    buffer[1] = (value & 0xFF00) >> 8;
    fwrite(buffer, 2, 1, file);
}

//Code inspired from the file in Course Documents in the Audio FIle Formats Test Tone Sample Code C file
void writeWavHeader(FILE *outputFile, WavHeader header , int bytes_per_sample, int numberSamples) {
    int dataChunkSize = 1 * numberSamples * bytes_per_sample;
    int formSize = 36 + dataChunkSize;
    short int frameSize = bytes_per_sample * header.num_channels;
    int bytesPerSecond =header.sample_rate * frameSize;
    fputs("RIFF", outputFile);
    fwriteIntLSB(formSize, outputFile);
    fputs("WAVE", outputFile);
    fputs("fmt ", outputFile);
    fwriteIntLSB(16, outputFile);
    fwriteShortLSB(1, outputFile);
    fwriteShortLSB((short)1, outputFile);
    fwriteIntLSB(header.sample_rate, outputFile);
    fwriteIntLSB(bytesPerSecond, outputFile);
    fwriteShortLSB(frameSize, outputFile);
    fwriteShortLSB(8 * bytes_per_sample, outputFile);
    fputs("data", outputFile);
    fwriteIntLSB(dataChunkSize, outputFile);
}

float shortToFloat(short value) {
   return value / 32768.0;
}

// Function to convolve two signals using the FFT
// Code inspired from Kimiya tutorial slides on CPSC501_F23_reverse_audio_time_convolution_FFT
void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

double* multiplyFrequencyData(double* freqData1, double* freqData2, int size) {
    double* result = (double*)malloc(size * sizeof(double));
    for (int i = 0; i < size; i += 2) {
        result[i] = freqData1[i] * freqData2[i] - freqData1[i+1] * freqData2[i+1]; // real part
        result[i + 1] = freqData1[i+1] * freqData2[i] + freqData1[i] * freqData2[i+1]; // imaginary part
    }
    return result;
}

double* convertToComplex(float* data, int dataSize, int arraySize) {
    double* complexData = (double*)malloc(arraySize * 2 * sizeof(double));

    for (int i = 0; i < arraySize; i++) {
        complexData[i] = 0.0;
    }

    for (int i = 0; i < dataSize; i++) {
        complexData[i * 2] = data[i]; // real part
    }
    return complexData;
}

float* convertToReal(double* complexData, int size) {
    float* result = (float*)malloc(size * sizeof(float));
    for (int i = 0; i < size; i += 2) {
        double real = complexData[i];
        double imag = complexData[i + 1];
        result[i / 2] = sqrt(real * real + imag * imag); // take the magnitude
    }
    return result;
}

// Function to write float data to a WAV file
void writeData(FILE *file, float data[], int size) {
    // Find the maximum absolute value
    float maxVal = 0.0f;
    for (int i = 0; i < size; ++i) {
        float absVal = fabs(data[i]);
        if (absVal > maxVal) {
            maxVal = absVal;
        }
    }

    // Normalize the data and write to file
    for (int i = 0; i < size; ++i) {
        // Normalize the data
        data[i] /= maxVal;
        // Convert to short and write to file
        short value = (short)(data[i] * 32768.0);
        fwrite(&value, sizeof(value), 1, file);

    }
}


/**
Read the tones, and call convolve on them
*/
//Code inspired from TA Ali Week 10 - Session 2 - Updated files
void readTone(char *sampleTone, char *impulseTone, char *output) {
    FILE *sampleFileStream = fopen(sampleTone, "rb");
    FILE *impulseFileStream = fopen(impulseTone, "rb");
    FILE *outputFileStream = fopen(output, "wb");

    WavHeader header_sample;
    WavHeader header_impulse;

    fread(&header_sample, sizeof(header_sample), 1, sampleFileStream);
    fread(&header_impulse, sizeof(header_impulse), 1, impulseFileStream);

    if (header_sample.subchunk1_size != 16) {
        int remainder = header_sample.subchunk1_size - 16;
        char randomVar[remainder];
        fread(randomVar, remainder, 1, sampleFileStream);
    }

    if (header_impulse.subchunk1_size != 16) {
        int remainder = header_impulse.subchunk1_size - 16;
        char randomVar[remainder];
        fread(randomVar, remainder, 1, impulseFileStream);
    }

    char subchunk2_id_sample[4];
    char subchunk2_id_impulse[4];
    int subchunk2_size_sample; 
    int subchunk2_size_impulse; 
    fread(&subchunk2_id_sample, sizeof(subchunk2_id_sample), 1, sampleFileStream);
    fread(&subchunk2_size_sample, sizeof(subchunk2_size_sample), 1, sampleFileStream);
    fread(&subchunk2_id_impulse, sizeof(subchunk2_id_impulse), 1, impulseFileStream);
    fread(&subchunk2_size_impulse, sizeof(subchunk2_size_impulse), 1, impulseFileStream);

    
    int bytesPerSample = header_sample.bits_per_sample / 8; // number of data points in the sample
    int num_samples = subchunk2_size_sample / bytesPerSample;
    int duration = num_samples / header_sample.sample_rate;
    int totalSamples = header_sample.sample_rate * duration;
 
    int bytesPerImpulse = header_impulse.bits_per_sample / 8; // number of data points in the sample
    int num_impulse = subchunk2_size_impulse / bytesPerImpulse;
    duration = num_impulse / header_impulse.sample_rate;
    int totalImpulses = header_impulse.sample_rate * duration;

    float *sampleData = (float *)malloc(totalSamples * sizeof(float));
    float *impulseData = (float *)malloc(totalImpulses * sizeof(float));

   for (int i = 0; i < totalSamples; ++i) { 
    short sampleValue;
    fread(&sampleValue, sizeof(sampleValue), 1, sampleFileStream);
    sampleData[i] = shortToFloat(sampleValue);
}

for (int i = 0; i < totalImpulses; ++i) { 
    short impulseValue;
    fread(&impulseValue, sizeof(impulseValue), 1, impulseFileStream);
    impulseData[i] = shortToFloat(impulseValue);
}

    // Convert the real data to complex data
    int maxSize = (int)pow(2, 24);
    double *complexSampleData = convertToComplex(sampleData, totalSamples, maxSize);
    double *complexImpulseData = convertToComplex(impulseData, totalImpulses, maxSize);

    // Close the input files
    fclose(sampleFileStream);
    fclose(impulseFileStream);

    // Convolve the data
    four1(complexSampleData - 1, totalSamples / 2, 1);
    four1(complexImpulseData - 1, totalImpulses / 2, 1);

    double *complexOutputData = multiplyFrequencyData(complexSampleData, complexImpulseData, maxSize);

    four1(complexOutputData - 1, maxSize / 2, -1);

    float *outputData = convertToReal(complexOutputData, maxSize);

    // Write the output WAV file
    int outputSize = totalSamples + totalImpulses - 1;
    writeWavHeader(outputFileStream, header_sample, bytesPerSample, outputSize);
    writeData(outputFileStream, outputData, outputSize);

    fclose(outputFileStream);

    // Clean up
    free(sampleData);
    free(impulseData);
    free(outputData);
    free(complexSampleData);
    free(complexImpulseData);
    free(complexOutputData);
}

// main line of execution
int main(int argc, char *argv[]) {
    char *sampleTone = NULL;
    char *impulseTone = NULL;
    char *output = NULL;

   
    if (argc == 4) {
        sampleTone = argv[1];
        impulseTone = argv[2];
        output = argv[3];
    } else {
        fprintf(stderr, "Usage:  %s sampleTone impulseTone\n", argv[0]);
        exit(-1);
    }
    readTone(sampleTone, impulseTone, output);
}
