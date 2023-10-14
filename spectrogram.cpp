#include <iostream>
#include <complex>
#include <vector>
#include <fftw3.h>
#include <sndfile.h>
#include <cstdlib>
#include <fstream>

// Define the values for M_PI, fftSize, and hopSize
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

const int fftSize = 1024; // Adjust the FFT size as needed
const int hopSize = 512; // Adjust the hop size as needed

// Function to create a spectrogram with time-frequency dimensions
std::vector<std::vector<double> > createSpectrogram(const char* inputWavFile, double& samplingFreq) {
    // Open the WAV file
    SF_INFO sfInfo;
    SNDFILE* sndfile = sf_open(inputWavFile, SFM_READ, &sfInfo);

    /*if (!sndfile) {
        std::cerr << "Error: Failed to open the input WAV file." << std::endl;
        exit(1);
    }*/

    samplingFreq = static_cast<double>(sfInfo.samplerate);

    // Parameters for the spectrogram
    int fftSize = 1024;  // Adjust the FFT size as needed
    int hopSize = 512;   // Adjust the hop size as needed

    // Create FFTW plans
    fftw_complex* fftInput = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * fftSize));
    fftw_complex* fftOutput = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * fftSize));
    fftw_plan fftPlan = fftw_plan_dft_1d(fftSize, fftInput, fftOutput, FFTW_FORWARD, FFTW_ESTIMATE);

    // Create a buffer for reading audio data
    std::vector<double> audioBuffer(fftSize, 0);

    // Create a vector to store the spectrogram
    std::vector<std::vector<double> > spectrogram;

    // Read and process the audio file
    while (sf_readf_double(sndfile, audioBuffer.data(), fftSize) > 0) {
        // Apply a window function to reduce spectral leakage (e.g., Hamming window)
        for (int i = 0; i < fftSize; i++) {
            audioBuffer[i] *= 0.54 - 0.46 * cos((float)(2 * M_PI * i / (fftSize - 1)));
        }

        // Perform FFT on the audio data
        for (int i = 0; i < fftSize; i++) {
            fftInput[i][0] = audioBuffer[i];
            fftInput[i][1] = 0.0;
        }

        fftw_execute(fftPlan);

        // Calculate magnitude of FFT and add to the spectrogram
        std::vector<double> magnitudes(fftSize / 2, 0);
        for (int i = 0; i < fftSize / 2; i++) {
            magnitudes[i] = std::abs(fftOutput[i][0]);
        }

        spectrogram.push_back(magnitudes);

        // Shift the input buffer for the next iteration
        for (int i = 0; i < fftSize - hopSize; i++) {
            audioBuffer[i] = audioBuffer[i + hopSize];
        }
    }

    // Clean up
    sf_close(sndfile);
    fftw_destroy_plan(fftPlan);
    fftw_free(fftInput);
    fftw_free(fftOutput);

    return spectrogram;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.wav output.txt" << std::endl;
        return 1;
    }

    const char* inputWavFile = argv[1];
    const char* outputTextFile = argv[2];

    double samplingFreq;
    // Create the spectrogram
    std::vector<std::vector<double> > spectrogram = createSpectrogram(inputWavFile, samplingFreq);

    // Save the spectrogram to a text file
    std::ofstream outputFile(outputTextFile);
    if (outputFile.is_open()) {
        // Write the spectrogram data along with time and frequency information
        for (size_t i = 0; i < spectrogram.size(); i++) {
            double timeInSeconds = static_cast<double>(i * hopSize) / samplingFreq;
            for (size_t j = 0; j < spectrogram[i].size(); j++) {
                double freqInHz = static_cast<double>(j) / fftSize * samplingFreq;
                outputFile << timeInSeconds << " " << freqInHz << " " << spectrogram[i][j] << std::endl;
            }
        }
        outputFile.close();
    } else {
        std::cerr << "Error: Failed to open the output file for writing." << std::endl;
        return 1;
    }

    std::cout << "Spectrogram saved to " << outputTextFile << std::endl;

    return 0;
}
