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

void printWavHeader(WavHeader header){
    printf("chunk_id: %.4s\n", header.chunk_id);
    printf("chunk_size: %d\n", header.chunk_size);
    printf("format: %.4s\n", header.format);
    printf("subchunk1_id: %.4s\n", header.subchunk1_id);
    printf("subchunk1_size: %d\n", header.subchunk1_size);
    printf("audio_format: %d\n", header.audio_format);
    printf("num_channels: %d\n", header.num_channels);
    printf("sample_rate: %d\n", header.sample_rate);
    printf("byte_rate: %d\n", header.byte_rate);
    printf("block_align: %d\n", header.block_align);
    printf("bits_per_sample: %d\n", header.bits_per_sample);
}

void fwriteIntLSB(int value, FILE *file) {
    unsigned char buffer[4];
    buffer[0] = (value & 0x000000FF);
    buffer[1] = (value & 0x0000FF00) >> 8;
    buffer[2] = (value & 0x00FF0000) >> 16;
    buffer[3] = (value & 0xFF000000) >> 24;
    fwrite(buffer, 4, 1, file);
}

void fwriteShortLSB(short value, FILE *file) {
    unsigned char buffer[2];
    buffer[0] = (value & 0x00FF);
    buffer[1] = (value & 0xFF00) >> 8;
    fwrite(buffer, 2, 1, file);
}


void writeWavHeader(FILE *outputFile, WavHeader header , int bytes_per_sample, int numberSamples) {
    int dataChunkSize = 1 * numberSamples * bytes_per_sample;
    int formSize = 36 + dataChunkSize;
    short int frameSize = 1 * bytes_per_sample;
    int bytesPerSecond = (int)ceil(header.sample_rate * frameSize);
    fputs("RIFF", outputFile);
    fwriteIntLSB(formSize, outputFile);
    fputs("WAVE", outputFile);
    fputs("fmt ", outputFile);
    fwriteIntLSB(16, outputFile);
    fwriteShortLSB(1, outputFile);
    fwriteShortLSB((short)1, outputFile);
    fwriteIntLSB((int)header.sample_rate, outputFile);
    fwriteIntLSB(bytesPerSecond, outputFile);
    fwriteShortLSB(frameSize, outputFile);
    fwriteShortLSB(bytes_per_sample, outputFile);
    fputs("data", outputFile);
    fwriteIntLSB(dataChunkSize, outputFile);
}

float shortToFloat(short value) {
   return value / 32768.0;
}

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
    FILE *textFile = fopen("output.txt", "w");
    for (int i = 0; i < size; ++i) {
        short value = (short)(data[i] * 32768.0);
        if (value > 1)
        value = 1;

        if (value < -1)
        value = -1;
        fwrite(&value, sizeof(value), 1, file);

        fprintf(textFile, "%hd\n", value);
    }
    fclose(textFile);
}


/**
Read the tones, and call convolve on them
*/
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

    
    int bytesPerSample = header_sample.bits_per_sample / 8; 
    int num_samples = subchunk2_size_sample / bytesPerSample;
    size_t num_size = subchunk2_size_sample; 

 
    int bytesPerImpulse = header_impulse.bits_per_sample / 8; 
    int num_impulse = subchunk2_size_impulse / bytesPerImpulse;
    size_t impulse_size = subchunk2_size_impulse; 

    float *sampleData = (float *)malloc(num_size * sizeof(float));
    float *impulseData = (float *)malloc(impulse_size * sizeof(float));

   for (int i = 0; i < num_samples; ++i) { 
    short sampleValue;
    fread(&sampleValue, sizeof(sampleValue), 1, sampleFileStream);
    sampleData[i] = sampleValue;
}

   for (int i = 0; i < num_impulse; ++i) { 
    short impulseValue;
    fread(&impulseValue, sizeof(impulseValue), 1, impulseFileStream);
    impulseData[i] = impulseValue;
}

    // Close the input files
    fclose(sampleFileStream);
    fclose(impulseFileStream);

    // Convolve the data
    int outputSize = num_samples + num_impulse - 1;
    float *outputData = (float *)malloc(outputSize * sizeof(float));
    convolve(sampleData, num_samples, impulseData, num_impulse, outputData, outputSize);

    // Write the output WAV file
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
