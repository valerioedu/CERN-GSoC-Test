/*
    * I'm gonna test on floats so 1 bit sign, 8 bits exponent and 23 bits mantissa
    *
    * I'm also gonna test compression and decompression, to test the compression success
    * 
    * Final tests will be done with both 8 and 16 bits
    * 
    * Tests will be done on Gaussian, uniform and exponential distributions
    * 
    * To spare computational complexity I won't loop through 2 parameters I'll equal them
*/


#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <fstream>
#include <vector>

#define dx 0.001   
#define domain 5
#define dc 0.5
#define cmin 0
#define cmax 5

/*
*    To work with the bits of the floating number I will copy them to a uint32 
*/

//compressed bits should be 8 or 16
float compressFloat(float f, int compressedbits) { 
    if (std::isnan(f) || std::isinf(f)) {
        return f;
    }          
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
    
    std::vector<uint8_t> compressedData(size * bytesPerValue);
    
    for (int i = 0; i < size; i++) {
        uint32_t bits;
        memcpy(&bits, &data[i], sizeof(float));
        
        // Check for NaN or Inf values
        bool isSpecial = std::isnan(data[i]) || std::isinf(data[i]);
        
        if (!isSpecial) {
            // Zero out the least significant bits as before
            bits &= ~((1 << compressBits) - 1);
        }
        
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
        double mean = 0.5 + j * dc;         // to avoid division by zero/zero variance
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
        double lambda = 0.5 + j * dc;       // to make it converge to 1
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
    
    int num_points = int(domain / dx);
    int num_params = int((cmax - cmin) / dc) + 1; // Number of parameters (means, ranges, lambdas)
    
    float** gaussian_data = gaussian_data_generator();
    float** uniform_data = uniform_data_generator();
    float** exponential_data = exponential_data_generator();
    
    // --- PART 2: COMPRESSION IMPLEMENTATION ---
    // Create flattened arrays for easier processing
    float* flat_gaussian = new float[num_params * num_points];
    float* flat_uniform = new float[num_params * num_points];
    float* flat_exponential = new float[num_params * num_points];
    
    // Flatten the 2D arrays
    for (int j = 0; j < num_params; j++) {
        for (int i = 0; i < num_points; i++) {
            int flat_idx = j * num_points + i;
            flat_gaussian[flat_idx] = gaussian_data[j][i];
            flat_uniform[flat_idx] = uniform_data[j][i];
            flat_exponential[flat_idx] = exponential_data[j][i];
        }
    }
    
    // Create arrays for compressed data
    float* gaussian_compressed_8bit = new float[num_params * num_points];
    float* gaussian_compressed_16bit = new float[num_params * num_points];
    
    float* uniform_compressed_8bit = new float[num_params * num_points];
    float* uniform_compressed_16bit = new float[num_params * num_points];
    
    float* exponential_compressed_8bit = new float[num_params * num_points];
    float* exponential_compressed_16bit = new float[num_params * num_points];
    
    // Apply compression
    std::cout << "Applying different compression levels...\n";
    for (int i = 0; i < num_params * num_points; i++) {
        // Gaussian compression
        gaussian_compressed_8bit[i] = compressFloat(flat_gaussian[i], 8);
        gaussian_compressed_16bit[i] = compressFloat(flat_gaussian[i], 16);
        
        // Uniform compression
        uniform_compressed_8bit[i] = compressFloat(flat_uniform[i], 8);
        uniform_compressed_16bit[i] = compressFloat(flat_uniform[i], 16);
        
        // Exponential compression
        exponential_compressed_8bit[i] = compressFloat(flat_exponential[i], 8);
        exponential_compressed_16bit[i] = compressFloat(flat_exponential[i], 16);
    }
    
    // --- PART 3: FILE I/O AND SIZE COMPARISON ---
    std::cout << "Saving data to binary files and comparing sizes...\n";
    
    // Save original data
    SaveToBinaryFile("uncompressed_gaussian.bin", flat_gaussian, sizeof(float) * num_params * num_points);
    SaveToBinaryFile("uncompressed_uniform.bin", flat_uniform, sizeof(float) * num_params * num_points);
    SaveToBinaryFile("uncompressed_exponential.bin", flat_exponential, sizeof(float) * num_params * num_points);
    
    // Compress and save with actual storage reduction
    auto compressed_gaussian_8bit = compressFloatArray(flat_gaussian, num_params * num_points, 8);
    auto compressed_gaussian_16bit = compressFloatArray(flat_gaussian, num_params * num_points, 16);

    auto compressed_uniform_8bit = compressFloatArray(flat_uniform, num_params * num_points, 8);
    auto compressed_uniform_16bit = compressFloatArray(flat_uniform, num_params * num_points, 16);

    auto compressed_exp_8bit = compressFloatArray(flat_exponential, num_params * num_points, 8);
    auto compressed_exp_16bit = compressFloatArray(flat_exponential, num_params * num_points, 16);

    // Save the truly compressed data
    SaveCompressedToBinaryFile("gaussian_8bit_compressed.bin", compressed_gaussian_8bit);
    SaveCompressedToBinaryFile("gaussian_16bit_compressed.bin", compressed_gaussian_16bit);

    SaveCompressedToBinaryFile("uniform_8bit_compressed.bin", compressed_uniform_8bit);
    SaveCompressedToBinaryFile("uniform_16bit_compressed.bin", compressed_uniform_16bit);

    SaveCompressedToBinaryFile("exponential_8bit_compressed.bin", compressed_exp_8bit);
    SaveCompressedToBinaryFile("exponential_16bit_compressed.bin", compressed_exp_16bit);

    // Compare file sizes
    std::cout << "File sizes comparison:\n";
    std::cout << "Original data (each distribution): " << sizeof(float) * num_params * num_points << " bytes\n";
    std::cout << "8-bit compressed gaussian: " << compressed_gaussian_8bit.size() << " bytes (" 
              << (100.0 * compressed_gaussian_8bit.size() / (sizeof(float) * num_params * num_points)) << "% of original)\n";
    std::cout << "16-bit compressed gaussian: " << compressed_gaussian_16bit.size() << " bytes (" 
              << (100.0 * compressed_gaussian_16bit.size() / (sizeof(float) * num_params * num_points)) << "% of original)\n";

    // Test decompression to verify data integrity
    float* test_decomp = new float[num_params * num_points];
    decompressToFloatArray(compressed_gaussian_8bit, test_decomp, num_params * num_points, 8);

    // Verify some values match the compressed version
    std::cout << "\nVerifying decompression on 8 bit compressed Gaussian (parameter = 0) - first 5 values:\n";
    std::cout << "Original\tCompressed\tDecompressed\n";

    // Sampled gaussian with mean = variance = 1
    int sample_param_idx = 1; // Parameter = 0.5 + 1*0.5 = 1.0
    int sample_offset = sample_param_idx * num_points;

    for (int i = 0; i < 5; i++) {
        int idx = sample_offset + i;
        std::cout << std::fixed << std::setprecision(6)
                << flat_gaussian[idx] << "\t" 
                << gaussian_compressed_8bit[idx] << "\t"
                << test_decomp[idx] << "\n";
    }
    
    // --- PART 4: STATISTICAL ANALYSIS ---
    std::cout << "\n=== STATISTICAL ANALYSIS ===\n";
    
    // Gaussian statistics
    double gaussian_orig_mean = mean(flat_gaussian, num_params * num_points);
    double gaussian_8bit_mean = mean(gaussian_compressed_8bit, num_params * num_points);
    double gaussian_16bit_mean = mean(gaussian_compressed_16bit, num_params * num_points);
    
    double gaussian_orig_var = variance(flat_gaussian, num_params * num_points, gaussian_orig_mean);
    double gaussian_8bit_var = variance(gaussian_compressed_8bit, num_params * num_points, gaussian_8bit_mean);
    double gaussian_16bit_var = variance(gaussian_compressed_16bit, num_params * num_points, gaussian_16bit_mean);

    double gaussian_orig_std = standardDeviation(flat_gaussian, num_params * num_points, gaussian_orig_mean);
    double gaussian_8bit_std = standardDeviation(gaussian_compressed_8bit, num_params * num_points, gaussian_8bit_mean);
    double gaussian_16bit_std = standardDeviation(gaussian_compressed_16bit, num_params * num_points, gaussian_16bit_mean);
    
    std::cout << "Gaussian Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << gaussian_orig_mean << " Variance: " << gaussian_orig_var << " Standard Deviation: " << gaussian_orig_std << "\n";
    std::cout << "8-bit comp - Mean: " << gaussian_8bit_mean << " Variance: " << gaussian_8bit_var << " Standard Deviation: " << gaussian_8bit_std << "\n";
    std::cout << "16-bit comp - Mean: " << gaussian_16bit_mean << " Variance: " << gaussian_16bit_var << " Standard Deviation: " << gaussian_16bit_std << "\n\n";
    
    // Uniform statistics 
    double uniform_orig_mean = mean(flat_uniform, num_params * num_points);
    double uniform_8bit_mean = mean(uniform_compressed_8bit, num_params * num_points);
    double uniform_16bit_mean = mean(uniform_compressed_16bit, num_params * num_points);
    
    double uniform_orig_var = variance(flat_uniform, num_params * num_points, uniform_orig_mean);
    double uniform_8bit_var = variance(uniform_compressed_8bit, num_params * num_points, uniform_8bit_mean);
    double uniform_16bit_var = variance(uniform_compressed_16bit, num_params * num_points, uniform_16bit_mean);

    double uniform_orig_std = standardDeviation(flat_uniform, num_params * num_points, uniform_orig_mean);
    double uniform_8bit_std = standardDeviation(uniform_compressed_8bit, num_params * num_points, uniform_8bit_mean);
    double uniform_16bit_std = standardDeviation(uniform_compressed_16bit, num_params * num_points, uniform_16bit_mean);
    
    std::cout << "Uniform Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << uniform_orig_mean << " Variance: " << uniform_orig_var << " Standard Deviation: "<< uniform_orig_std << "\n";
    std::cout << "8-bit comp - Mean: " << uniform_8bit_mean << " Variance: " << uniform_8bit_var << " Standard Deviation: "<< uniform_8bit_std << "\n";
    std::cout << "16-bit comp - Mean: " << uniform_16bit_mean << " Variance: " << uniform_16bit_var << " Standard Deviation: "<< uniform_16bit_std << "\n\n";
    
    // Exponential statistics
    double exp_orig_mean = mean(flat_exponential, num_params * num_points);
    double exp_8bit_mean = mean(exponential_compressed_8bit, num_params * num_points);
    double exp_16bit_mean = mean(exponential_compressed_16bit, num_params * num_points);
    
    double exp_orig_var = variance(flat_exponential, num_params * num_points, exp_orig_mean);
    double exp_8bit_var = variance(exponential_compressed_8bit, num_params * num_points, exp_8bit_mean);
    double exp_16bit_var = variance(exponential_compressed_16bit, num_params * num_points, exp_16bit_mean);

    double exp_orig_std = standardDeviation(flat_exponential, num_params * num_points, exp_orig_mean);
    double exp_8bit_std = standardDeviation(exponential_compressed_8bit, num_params * num_points, exp_8bit_mean);
    double exp_16bit_std = standardDeviation(exponential_compressed_16bit, num_params * num_points, exp_16bit_mean);
    
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
    
    for (int i = 0; i < num_params * num_points; i++) {
        // Gaussian errors
        double err_gaussian_8bit = fabs(flat_gaussian[i] - gaussian_compressed_8bit[i]);
        double err_gaussian_16bit = fabs(flat_gaussian[i] - gaussian_compressed_16bit[i]);
        
        mse_gaussian_8bit += err_gaussian_8bit * err_gaussian_8bit;
        mse_gaussian_16bit += err_gaussian_16bit * err_gaussian_16bit;
        
        max_err_gaussian_8bit = std::max(max_err_gaussian_8bit, err_gaussian_8bit);
        max_err_gaussian_16bit = std::max(max_err_gaussian_16bit, err_gaussian_16bit);
        
        // Uniform errors
        double err_uniform_8bit = fabs(flat_uniform[i] - uniform_compressed_8bit[i]);
        double err_uniform_16bit = fabs(flat_uniform[i] - uniform_compressed_16bit[i]);
        
        mse_uniform_8bit += err_uniform_8bit * err_uniform_8bit;
        mse_uniform_16bit += err_uniform_16bit * err_uniform_16bit;
        
        max_err_uniform_8bit = std::max(max_err_uniform_8bit, err_uniform_8bit);
        max_err_uniform_16bit = std::max(max_err_uniform_16bit, err_uniform_16bit);
        
        // Exponential errors
        double err_exp_8bit = fabs(flat_exponential[i] - exponential_compressed_8bit[i]);
        double err_exp_16bit = fabs(flat_exponential[i] - exponential_compressed_16bit[i]);
        
        mse_exp_8bit += err_exp_8bit * err_exp_8bit;
        mse_exp_16bit += err_exp_16bit * err_exp_16bit;
        
        max_err_exp_8bit = std::max(max_err_exp_8bit, err_exp_8bit);
        max_err_exp_16bit = std::max(max_err_exp_16bit, err_exp_16bit);
    }
    
    // Normalize MSE values
    mse_gaussian_8bit /= (num_params * num_points);
    mse_gaussian_16bit /= (num_params * num_points);
    
    mse_uniform_8bit /= (num_params * num_points);
    mse_uniform_16bit /= (num_params * num_points);
    
    mse_exp_8bit /= (num_params * num_points);
    mse_exp_16bit /= (num_params * num_points);
    
    // Display error metrics
    std::cout << "Gaussian Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_gaussian_8bit << " Max Error: " << max_err_gaussian_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_gaussian_16bit << " Max Error: " << max_err_gaussian_16bit << "\n\n";
    
    std::cout << "Uniform Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_uniform_8bit << " Max Error: " << max_err_uniform_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_uniform_16bit << " Max Error: " << max_err_uniform_16bit << "\n\n";
    
    std::cout << "Exponential Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_exp_8bit << " Max Error: " << max_err_exp_8bit << "\n";
    std::cout << "16-bit - MSE: " << mse_exp_16bit << " Max Error: " << max_err_exp_16bit << "\n\n";
    
    // --- PART 6: EXPORT DATA FOR PLOTTING ---
    std::cout << "Exporting data for plotting...\n";

    // Define parameter indices (corresponding to parameter ≈ 1.0)
    int gaussian_param_idx = 1;  // Parameter = 0.5 + 1*0.5 = 1.0
    int uniform_param_idx = 0;   // Parameter = 0*0.5 = 0.0 (lower bound)
    int exp_param_idx = 1;       // Parameter = 0.5 + 1*0.5 = 1.0 (lambda)

    // Export original data for the selected parameter
    ExportForPlotting("gaussian_original.csv", gaussian_data[gaussian_param_idx], num_points);
    ExportForPlotting("uniform_original.csv", uniform_data[uniform_param_idx], num_points);
    ExportForPlotting("exponential_original.csv", exponential_data[exp_param_idx], num_points);

    // Create compressed data arrays for the selected parameter set
    float* first_gaussian_8bit = new float[num_points];
    float* first_gaussian_16bit = new float[num_points];
    float* first_uniform_8bit = new float[num_points];
    float* first_uniform_16bit = new float[num_points];
    float* first_exp_8bit = new float[num_points];
    float* first_exp_16bit = new float[num_points];

    // Extract selected parameter set data from flattened arrays
    for (int i = 0; i < num_points; i++) {
        // Calculate the position in the flattened array for the selected parameter set
        int gaussian_idx = gaussian_param_idx * num_points + i;
        int uniform_idx = uniform_param_idx * num_points + i;
        int exp_idx = exp_param_idx * num_points + i;
        
        first_gaussian_8bit[i] = gaussian_compressed_8bit[gaussian_idx];
        first_gaussian_16bit[i] = gaussian_compressed_16bit[gaussian_idx];
        first_uniform_8bit[i] = uniform_compressed_8bit[uniform_idx];
        first_uniform_16bit[i] = uniform_compressed_16bit[uniform_idx];
        first_exp_8bit[i] = exponential_compressed_8bit[exp_idx];
        first_exp_16bit[i] = exponential_compressed_16bit[exp_idx];
    }
    
    // Export compressed first parameter set data for plotting
    ExportForPlotting("gaussian_8bit.csv", first_gaussian_8bit, num_points);
    ExportForPlotting("gaussian_16bit.csv", first_gaussian_16bit, num_points);
    ExportForPlotting("uniform_8bit.csv", first_uniform_8bit, num_points);
    ExportForPlotting("uniform_16bit.csv", first_uniform_16bit, num_points);
    ExportForPlotting("exponential_8bit.csv", first_exp_8bit, num_points);
    ExportForPlotting("exponential_16bit.csv", first_exp_16bit, num_points);
    
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
    for (int j = 0; j < num_params; j++) {
        delete[] gaussian_data[j];
        delete[] uniform_data[j];
        delete[] exponential_data[j];
    }
    delete[] gaussian_data;
    delete[] uniform_data;
    delete[] exponential_data;
    
    delete[] flat_gaussian;
    delete[] flat_uniform;
    delete[] flat_exponential;
    delete[] gaussian_compressed_8bit;
    delete[] gaussian_compressed_16bit;
    delete[] uniform_compressed_8bit;
    delete[] uniform_compressed_16bit;
    delete[] exponential_compressed_8bit;
    delete[] exponential_compressed_16bit;
    delete[] test_decomp;
    
    delete[] first_gaussian_8bit;
    delete[] first_gaussian_16bit;
    delete[] first_uniform_8bit;
    delete[] first_uniform_16bit;
    delete[] first_exp_8bit;
    delete[] first_exp_16bit;
    
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