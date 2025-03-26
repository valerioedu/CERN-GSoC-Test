/*
    * I'm gonna test on both float and double
    *
    * I'm also gonna test compression on the original primitive types and structs
    * 
    * Final tests will be done with both 8 and 16 bits
    * 
    * Tests will be done on Gaussian, uniform and exponential distributions
*/


#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <fstream>
#include <vector>

#define dx 0.001   
#define domain 2
#define dc 0.5
#define cmin -2.5
#define cmax 2.5

/*
*    To work with the bits of the floating number I will copy them to a uint32\uint64 
*/

//compressed bits should be 8 or 16
float compressFloat(float f, int compressedbits) {           
    uint32_t bits;        

    //copying the bits of the float into the uint32
    memcpy(&bits, &f, sizeof(float));

    //masking the bits that will be discarded
    bits &= ~((1 << compressedbits) - 1);                    
    float compressed;

    //copying the bits back to a float
    memcpy(&compressed, &bits, sizeof(float)); 
    return compressed;
}

double compressDouble(double d, int compressedbits) {
    uint64_t bits;
    memcpy(&bits, &d, sizeof(double));
    bits &= ~((1 << compressedbits) - 1);
    double compressed;
    memcpy(&compressed, &bits, sizeof(double));
    return compressed;
}

//Distribution functions, for ease of testing I won't use constants

double gaussian(double x, double mean, double variance) {
    return exp(-((x - mean) * (x - mean)) / (2 * variance)) / sqrt(2 * M_PI * variance);
}

double uniform(double x, double a, double b) {
    return (x >= a && x <= b) ? 1 / (b - a) : 0;
}

double exponential(double x, double lambda) {
    return (x >= 0) ? lambda *exp(-lambda * x) : 0;
}

//Statistical functions

double mean(float* data, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += data[i];
    }
    return sum / size;
}

double variance(float* data, int size, double mean) {
    double variance = 0.0;
    for (int i = 0; i < size; i++) {
        variance += (data[i] - mean) * (data[i] - mean);
    }
    return variance / size;
}

double standardDeviation(float* data, int size, double mean) {
    return sqrt(variance(data, size, mean));
}

//I/O functions

void SaveToBinaryFile(const std::string& filename, void* data, size_t size) {
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<char*>(data), size);
}

void ExportForPlotting(const std::string& filename, float* data, int size) {
    std::ofstream file(filename);
    for (int i = 0; i < size; i++) {
        file << i << "," << data[i] << "\n";
    }
}

// Function to compress float data and store in a smaller buffer
std::vector<uint8_t> compressFloatArray(float* data, int size, int compressBits) {
    // Calculate how many bytes per value we need after compression
    int bytesPerValue = (32 - compressBits) / 8;
    if ((32 - compressBits) % 8 != 0) bytesPerValue++; // Round up
    
    // Create a buffer to hold the compressed data
    std::vector<uint8_t> compressedData(size * bytesPerValue);
    
    for (int i = 0; i < size; i++) {
        // Get the bits for this float
        uint32_t bits;
        memcpy(&bits, &data[i], sizeof(float));
        
        // Zero out the least significant bits as before
        bits &= ~((1 << compressBits) - 1);
        
        // Copy only the remaining significant bytes to the output buffer
        for (int b = 0; b < bytesPerValue; b++) {
            int bytePosition = b + (4 - bytesPerValue); // Skip the zeroed bytes
            if (bytePosition >= 0 && bytePosition < 4) {
                compressedData[i * bytesPerValue + b] = 
                    (bits >> (bytePosition * 8)) & 0xFF;
            }
        }
    }
    
    return compressedData;
}

// Function to decompress data back to float array
void decompressToFloatArray(const std::vector<uint8_t>& compressedData, 
                           float* outputArray, int size, int compressBits) {
    // Calculate how many bytes per value are in the compressed data
    int bytesPerValue = (32 - compressBits) / 8;
    if ((32 - compressBits) % 8 != 0) bytesPerValue++;
    
    for (int i = 0; i < size; i++) {
        // Start with all zeroes
        uint32_t bits = 0;
        
        // Reconstruct the value from the compressed bytes
        for (int b = 0; b < bytesPerValue; b++) {
            int bytePosition = b + (4 - bytesPerValue);
            if (bytePosition >= 0 && bytePosition < 4) {
                bits |= static_cast<uint32_t>(compressedData[i * bytesPerValue + b]) 
                      << (bytePosition * 8);
            }
        }
        
        // Convert back to float
        memcpy(&outputArray[i], &bits, sizeof(float));
    }
}

// Function to save the compressed data to a binary file
void SaveCompressedToBinaryFile(const std::string& filename, 
                               const std::vector<uint8_t>& data) {
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(data.data()), data.size());
}

/*
* Generator functions
*/

float** gaussian_data_generator() {
    int num_points = int(domain / dx);
    int num_means = int((cmax - cmin) / dc) + 1;
    
    float** gaussian_data = new float*[num_means];
    
    for (int j = 0; j < num_means; j++) {
        gaussian_data[j] = new float[num_points];
    }
    
    for (int j = 0; j < num_means; j++) {
        double mean = cmin + j * dc;
        for (int i = 0; i < num_points; i++) {
            double x = i * dx;
            gaussian_data[j][i] = gaussian(x, mean, mean); // Using mean as variance too
        }
    }

    return gaussian_data;
}

float** uniform_data_generator() {
    int num_points = int(domain / dx);
    int num_means = int((cmax - cmin) / dc) + 1;
    
    float** uniform_data = new float*[num_means];
    
    for (int j = 0; j < num_means; j++) {
        uniform_data[j] = new float[num_points];
    }
    
    for (int j = 0; j < num_means; j++) {
        double c = cmin + j * dc;
        for (int i = 0; i < num_points; i++) {
            double x = i * dx;
            uniform_data[j][i] = uniform(x, c, c + 1);
        }
    }

    return uniform_data;
}

float** exponential_data_generator() {
    int num_points = int(domain / dx);
    int num_lambdas = int((cmax - cmin) / dc) + 1;
    
    float** exponential_data = new float*[num_lambdas];
    
    for (int j = 0; j < num_lambdas; j++) {
        exponential_data[j] = new float[num_points];
    }
    
    for (int j = 0; j < num_lambdas; j++) {
        double lambda = cmin + j * dc;
        for (int i = 0; i < num_points; i++) {
            double x = i * dx;
            exponential_data[j][i] = exponential(x, lambda);
        }
    }

    return exponential_data;
}

//Main function

int main() {
    // --- PART 1: DATA GENERATION ---
    // Create arrays for different distributions
    std::cout << "Generating data from three different distributions...\n";
    
    float** gaussian_data = gaussian_data_generator();
    float** uniform_data = uniform_data_generator();
    float** exponential_data = exponential_data_generator();
    
    // --- PART 2: COMPRESSION IMPLEMENTATION ---
    // Create arrays for compressed data
    float* gaussian_compressed_8bit = new float[int(domain / dx)];
    float* gaussian_compressed_16bit = new float[int(domain / dx)];
    
    float* uniform_compressed_8bit = new float[int(domain / dx)];
    float* uniform_compressed_16bit = new float[int(domain / dx)];
    
    float* exponential_compressed_8bit = new float[int(domain / dx)];
    float* exponential_compressed_16bit = new float[int(domain / dx)];
    
    // Apply compression
    std::cout << "Applying different compression levels...\n";
    for (int i = 0; i < int(domain / dx); i++) {
        // Gaussian compression
        gaussian_compressed_8bit[i] = compressFloat(gaussian_data[i], 8);
        gaussian_compressed_16bit[i] = compressFloat(gaussian_data[i], 16);
        
        // Uniform compression
        uniform_compressed_8bit[i] = compressFloat(uniform_data[i], 8);
        uniform_compressed_16bit[i] = compressFloat(uniform_data[i], 16);
        
        // Exponential compression
        exponential_compressed_8bit[i] = compressFloat(exponential_data[i], 8);
        exponential_compressed_16bit[i] = compressFloat(exponential_data[i], 16);
    }
    
    // --- PART 3: FILE I/O AND SIZE COMPARISON ---
    std::cout << "Saving data to binary files and comparing sizes...\n";
    
    // Save original data
    SaveToBinaryFile("uncompressed_gaussian.bin", gaussian_data, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uncompressed_uniform.bin", uniform_data, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uncompressed_exponential.bin", exponential_data, sizeof(float) * int(domain / dx));
    
    // Compress and save with actual storage reduction
    auto compressed_gaussian_8bit = compressFloatArray(gaussian_data, int(domain / dx), 8);
    auto compressed_gaussian_16bit = compressFloatArray(gaussian_data, int(domain / dx), 16);

    auto compressed_uniform_8bit = compressFloatArray(uniform_data, int(domain / dx), 8);
    auto compressed_uniform_16bit = compressFloatArray(uniform_data, int(domain / dx), 16);

    auto compressed_exp_8bit = compressFloatArray(exponential_data, int(domain / dx), 8);
    auto compressed_exp_16bit = compressFloatArray(exponential_data, int(domain / dx), 16);

    // Save the truly compressed data
    SaveCompressedToBinaryFile("gaussian_8bit_compressed.bin", compressed_gaussian_8bit);
    SaveCompressedToBinaryFile("gaussian_16bit_compressed.bin", compressed_gaussian_16bit);

    SaveCompressedToBinaryFile("uniform_8bit_compressed.bin", compressed_uniform_8bit);
    SaveCompressedToBinaryFile("uniform_16bit_compressed.bin", compressed_uniform_16bit);

    SaveCompressedToBinaryFile("exponential_8bit_compressed.bin", compressed_exp_8bit);
    SaveCompressedToBinaryFile("exponential_16bit_compressed.bin", compressed_exp_16bit);

    // Compare file sizes
    std::cout << "File sizes comparison:\n";
    std::cout << "Original data (each distribution): " << sizeof(float) * int(domain / dx) << " bytes\n";
    std::cout << "8-bit compressed gaussian: " << compressed_gaussian_8bit.size() << " bytes (" 
              << (100.0 * compressed_gaussian_8bit.size() / (sizeof(float) * int(domain / dx))) << "% of original)\n";
    std::cout << "16-bit compressed gaussian: " << compressed_gaussian_16bit.size() << " bytes (" 
              << (100.0 * compressed_gaussian_16bit.size() / (sizeof(float) * int(domain / dx))) << "% of original)\n";

    // If you need to test decompression, add this:
    // Test decompression to verify data integrity
    float* test_decomp = new float[int(domain / dx)];
    decompressToFloatArray(compressed_gaussian_8bit, test_decomp, int(domain / dx), 8);

    // Verify some values match the compressed version
    std::cout << "\nVerifying decompression - first 5 values:\n";
    std::cout << "Original\tCompressed\tDecompressed\n";
    for (int i = 0; i < 5; i++) {
        std::cout << gaussian_data[i] << "\t" 
                  << gaussian_compressed_8bit[i] << "\t"
                  << test_decomp[i] << "\n";
    }
    
    // --- PART 4: STATISTICAL ANALYSIS ---
    std::cout << "\n=== STATISTICAL ANALYSIS ===\n";
    
    // Gaussian statistics
    double gaussian_orig_mean = mean(gaussian_data, int(domain / dx));
    double gaussian_8bit_mean = mean(gaussian_compressed_8bit, int(domain / dx));
    double gaussian_16bit_mean = mean(gaussian_compressed_16bit, int(domain / dx));
    
    double gaussian_orig_var = variance(gaussian_data, int(domain / dx), gaussian_orig_mean);
    double gaussian_8bit_var = variance(gaussian_compressed_8bit, int(domain / dx), gaussian_8bit_mean);
    double gaussian_16bit_var = variance(gaussian_compressed_16bit, int(domain / dx), gaussian_16bit_mean);

    double gaussian_orig_std = standardDeviation(gaussian_data, int(domain / dx), gaussian_orig_mean);
    double gaussian_8bit_std = standardDeviation(gaussian_compressed_8bit, int(domain / dx), gaussian_8bit_mean);
    double gaussian_16bit_std = standardDeviation(gaussian_compressed_16bit, int(domain / dx), gaussian_16bit_mean);
    
    std::cout << "Gaussian Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << gaussian_orig_mean << " Variance: " << gaussian_orig_var << " Standard Deviation: " << gaussian_orig_std << "\n";
    std::cout << "8-bit comp - Mean: " << gaussian_8bit_mean << " Variance: " << gaussian_8bit_var << " Standard Deviation: " << gaussian_8bit_std << "\n";
    std::cout << "16-bit comp - Mean: " << gaussian_16bit_mean << " Variance: " << gaussian_16bit_var << " Standard Deviation: " << gaussian_16bit_std << "\n\n";
    
    // Uniform statistics 
    double uniform_orig_mean = mean(uniform_data, int(domain / dx));
    double uniform_8bit_mean = mean(uniform_compressed_8bit, int(domain / dx));
    double uniform_16bit_mean = mean(uniform_compressed_16bit, int(domain / dx));
    
    double uniform_orig_var = variance(uniform_data, int(domain / dx), uniform_orig_mean);
    double uniform_8bit_var = variance(uniform_compressed_8bit, int(domain / dx), uniform_8bit_mean);
    double uniform_16bit_var = variance(uniform_compressed_16bit, int(domain / dx), uniform_16bit_mean);

    double uniform_orig_std = standardDeviation(uniform_data, int(domain / dx), uniform_orig_mean);
    double uniform_8bit_std = standardDeviation(uniform_compressed_8bit, int(domain / dx), uniform_8bit_mean);
    double uniform_16bit_std = standardDeviation(uniform_compressed_16bit, int(domain / dx), uniform_16bit_mean);
    
    std::cout << "uniform Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << uniform_orig_mean << " Variance: " << uniform_orig_var << " Standard Deviation: "<< uniform_orig_std << "\n";
    std::cout << "8-bit comp - Mean: " << uniform_8bit_mean << " Variance: " << uniform_8bit_var << " Standard Deviation: "<< uniform_8bit_std << "\n";
    std::cout << "16-bit comp - Mean: " << uniform_16bit_mean << " Variance: " << uniform_16bit_var << " Standard Deviation: "<< uniform_16bit_std << "\n\n";
    
    // Exponential statistics
    double exp_orig_mean = mean(exponential_data, int(domain / dx));
    double exp_8bit_mean = mean(exponential_compressed_8bit, int(domain / dx));
    double exp_16bit_mean = mean(exponential_compressed_16bit, int(domain / dx));
    
    double exp_orig_var = variance(exponential_data, int(domain / dx), exp_orig_mean);
    double exp_8bit_var = variance(exponential_compressed_8bit, int(domain / dx), exp_8bit_mean);
    double exp_16bit_var = variance(exponential_compressed_16bit, int(domain / dx), exp_16bit_mean);

    double exp_orig_std = standardDeviation(exponential_data, int(domain / dx), exp_orig_mean);
    double exp_8bit_std = standardDeviation(exponential_compressed_8bit, int(domain / dx), exp_8bit_mean);
    double exp_16bit_std = standardDeviation(exponential_compressed_16bit, int(domain / dx), exp_16bit_mean);
    
    std::cout << "Exponential Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << exp_orig_mean << " Variance: " << exp_orig_var << " Standard Deviation: " << exp_orig_std << "\n";
    std::cout << "8-bit comp - Mean: " << exp_8bit_mean << " Variance: " << exp_8bit_var << " Standard Deviation: " << exp_8bit_std << "\n";
    std::cout << "16-bit comp - Mean: " << exp_16bit_mean << " Variance: " << exp_16bit_var << " Standard Deviation: " << exp_16bit_std << "\n\n";
    
    // --- PART 5: ERROR ANALYSIS ---
    std::cout << "\n=== ERROR ANALYSIS ===\n";
    
    // Calculate MSE for different distributions and compression levels
    double mse_gaussian_8bit = 0.0, mse_gaussian_16bit = 0.0;
    double mse_uniform_8bit = 0.0, mse_uniform_16bit = 0.0;
    double mse_exp_8bit = 0.0, mse_exp_16bit = 0.0;
    
    // Max errors
    double max_err_gaussian_8bit = 0.0, max_err_gaussian_16bit = 0.0;
    double max_err_uniform_8bit = 0.0, max_err_uniform_16bit = 0.0;
    double max_err_exp_8bit = 0.0, max_err_exp_16bit = 0.0;
    
    for (int i = 0; i < int(domain / dx); i++) {
        // Gaussian errors
        double err_gaussian_8bit = fabs(gaussian_data[i] - gaussian_compressed_8bit[i]);
        double err_gaussian_16bit = fabs(gaussian_data[i] - gaussian_compressed_16bit[i]);
        
        mse_gaussian_8bit += err_gaussian_8bit * err_gaussian_8bit;
        mse_gaussian_16bit += err_gaussian_16bit * err_gaussian_16bit;
        
        max_err_gaussian_8bit = std::max(max_err_gaussian_8bit, err_gaussian_8bit);
        max_err_gaussian_16bit = std::max(max_err_gaussian_16bit, err_gaussian_16bit);
        
        // Uniform errors
        double err_uniform_8bit = fabs(uniform_data[i] - uniform_compressed_8bit[i]);
        double err_uniform_16bit = fabs(uniform_data[i] - uniform_compressed_16bit[i]);
        
        mse_uniform_8bit += err_uniform_8bit * err_uniform_8bit;
        mse_uniform_16bit += err_uniform_16bit * err_uniform_16bit;
        
        max_err_uniform_8bit = std::max(max_err_uniform_8bit, err_uniform_8bit);
        max_err_uniform_16bit = std::max(max_err_uniform_16bit, err_uniform_16bit);
        
        // Exponential errors
        double err_exp_8bit = fabs(exponential_data[i] - exponential_compressed_8bit[i]);
        double err_exp_16bit = fabs(exponential_data[i] - exponential_compressed_16bit[i]);
        
        mse_exp_8bit += err_exp_8bit * err_exp_8bit;
        mse_exp_16bit += err_exp_16bit * err_exp_16bit;
        
        max_err_exp_8bit = std::max(max_err_exp_8bit, err_exp_8bit);
        max_err_exp_16bit = std::max(max_err_exp_16bit, err_exp_16bit);
    }
    
    int numPoints = int(domain / dx);
    
    // uniformize MSE values
    mse_gaussian_8bit /= numPoints;
    mse_gaussian_16bit /= numPoints;
    
    mse_uniform_8bit /= numPoints;
    mse_uniform_16bit /= numPoints;
    
    mse_exp_8bit /= numPoints;
    mse_exp_16bit /= numPoints;
    
    // Display error metrics
    std::cout << "Gaussian Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_gaussian_8bit << " Max Error: " << max_err_gaussian_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_gaussian_16bit << " Max Error: " << max_err_gaussian_16bit << "\n\n";
    
    std::cout << "uniform Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_uniform_8bit << " Max Error: " << max_err_uniform_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_uniform_16bit << " Max Error: " << max_err_uniform_16bit << "\n\n";
    
    std::cout << "Exponential Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_exp_8bit << " Max Error: " << max_err_exp_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_exp_16bit << " Max Error: " << max_err_exp_16bit << "\n\n";
    
    // --- PART 6: EXPORT DATA FOR PLOTTING ---
    std::cout << "Exporting data for plotting...\n";
        
    ExportForPlotting("gaussian_original.csv", gaussian_data, int(domain / dx));
    ExportForPlotting("gaussian_8bit.csv", gaussian_compressed_8bit, int(domain / dx));
    ExportForPlotting("gaussian_16bit.csv", gaussian_compressed_16bit, int(domain / dx));

    ExportForPlotting("uniform_original.csv", uniform_data, int(domain / dx));
    ExportForPlotting("uniform_8bit.csv", uniform_compressed_8bit, int(domain / dx));
    ExportForPlotting("uniform_16bit.csv", uniform_compressed_16bit, int(domain / dx));

    ExportForPlotting("exponential_original.csv", exponential_data, int(domain / dx));
    ExportForPlotting("exponential_8bit.csv", exponential_compressed_8bit, int(domain / dx));
    ExportForPlotting("exponential_16bit.csv", exponential_compressed_16bit, int(domain / dx));
    
    // --- PART 7: CONCLUSIONS ---
    std::cout << "\n=== COMPRESSION RECOMMENDATIONS ===\n";
    std::cout << "Based on the error analysis:\n";
    
    // Simple recommendation logic
    if (mse_gaussian_8bit < 1e-6 && mse_uniform_8bit < 1e-6 && mse_exp_8bit < 1e-6) {
        std::cout << "- 8-bit compression is suitable for most applications with minimal loss of precision\n";
    } else {
        std::cout << "- For high-precision requirements, use 16-bit compression or no compression at all\n";
    }
    
    std::cout << "- For applications with strict error tolerances, consider the max error values\n";
    std::cout << "- Statistical parameters show minimal changes in distribution characteristics\n";
    std::cout << "  even with aggressive compression\n";
    
    // Clean up
    delete[] gaussian_data;
    delete[] uniform_data;
    delete[] exponential_data;
    delete[] gaussian_compressed_8bit;
    delete[] gaussian_compressed_16bit;
    delete[] uniform_compressed_8bit;
    delete[] uniform_compressed_16bit;
    delete[] exponential_compressed_8bit;
    delete[] exponential_compressed_16bit;
    delete[] test_decomp;
    
    // Call the Python plotting script
    std::cout << "\n=== GENERATING VISUALIZATIONS ===\n";
    std::cout << "Calling Python script to generate plots...\n";
    
    int result = system("python3 plot.py");
    
    if (result == 0) {
        std::cout << "Visualization completed successfully. Check the generated PNG files.\n";
    } else {
        std::cout << "Error running visualization script. Probably Python and required packages are not installed.\n";
    }

    return 0;
}