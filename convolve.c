#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void writeFloatArrayToFile(const char *filename, float data[], int size) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; ++i) {
        fprintf(file, "%f\n", data[i]);
    }

    fclose(file);
}

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
    short int frameSize = 1 * bytes_per_sample;
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

//Slides taken from sources iwth the CPSC 501 Course Documents of TimeDomainConvolution
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n,m;

    for (n=0; n < P; n++)
    {
        y[n] = 0.0;
    }

    for (n=0; n<N; n++){
        for (m=0; m<M; m++){
            y[n+m] += x[n] * h[m];
        }
    }
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

    // Close the input files
    fclose(sampleFileStream);
    fclose(impulseFileStream);

    // Convolve the data
    int outputSize = totalSamples + totalImpulses - 1;
    float *outputData = (float *)malloc(outputSize * sizeof(float));
    convolve(sampleData, totalSamples, impulseData, totalImpulses, outputData, outputSize);

    // Write the output WAV file
    writeFloatArrayToFile("data.txt", outputData, outputSize);
    writeWavHeader(outputFileStream, header_sample, bytesPerSample, outputSize);
    writeData(outputFileStream, outputData, outputSize);

    fclose(outputFileStream);

    // Clean up
    free(sampleData);
    free(impulseData);
    free(outputData);
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
