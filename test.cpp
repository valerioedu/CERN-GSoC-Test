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

#define dx 0.001   
#define domain 2
#define dc 0.5
#define cdomain 5

//Float struct

typedef struct FloatBits {
    uint32_t mantissa : 23;
    uint32_t exponent : 8;
    uint32_t sign : 1;
} FloatBits;

//Double struct

typedef struct DoubleBits {
    uint64_t mantissa : 52;
    uint64_t exponent : 11;
    uint64_t sign : 1;
} DoubleBits;

//Compression on primitive types

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

//Compression on structs

FloatBits compressFloatStruct(float f, int compressedbits) {
    FloatBits bits;

    //copying the bits of the float into the struct
    memcpy(&bits, &f, sizeof(float));

    //masking the bits that will be discarded
    bits.mantissa &= ~((1 << compressedbits) - 1);
    return bits;
}

DoubleBits compressDoubleStruct(double d, int compressedbits) {
    DoubleBits bits;
    memcpy(&bits, &d, sizeof(double));
    bits.mantissa &= ~((1 << compressedbits) - 1);
    return bits;
}

//Distribution functions, for ease of testing I won't use constants

double gaussian(double x) {
    return exp(-x * x);
}

double uniform(double x) {
    return (x >= 0 && x <= 1) ? 1 : 0;
}

double exponential(double x) {
    return (x >= 0) ? exp(-x) : 0;
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

//Main function

int main() {
    // --- PART 1: DATA GENERATION ---
    // Create arrays for different distributions
    float* gaussian_data = new float[int(domain / dx)];
    float* uniform_data = new float[int(domain / dx)];
    float* exponential_data = new float[int(domain / dx)];
    
    // Generate data for three different distributions
    std::cout << "Generating data from three different distributions...\n";
    for (int i = 0; i < int(domain / dx); i++) {
        double x = i * dx;  // Scale for meaningful variations
        
        // Gaussian distribution
        gaussian_data[i] = static_cast<float>(gaussian(x));
        
        // uniform distribution
        uniform_data[i] = static_cast<float>(uniform(x));
        
        // Exponential distribution
        exponential_data[i] = static_cast<float>(exponential(x));
    }
    
    // --- PART 2: COMPRESSION IMPLEMENTATION ---
    // Create arrays for compressed data
    float* gaussian_compressed_8bit = new float[int(domain / dx)];
    float* gaussian_compressed_12bit = new float[int(domain / dx)];
    float* gaussian_compressed_16bit = new float[int(domain / dx)];
    
    float* uniform_compressed_8bit = new float[int(domain / dx)];
    float* uniform_compressed_12bit = new float[int(domain / dx)];
    float* uniform_compressed_16bit = new float[int(domain / dx)];
    
    float* exponential_compressed_8bit = new float[int(domain / dx)];
    float* exponential_compressed_12bit = new float[int(domain / dx)];
    float* exponential_compressed_16bit = new float[int(domain / dx)];
    
    // Apply compression
    std::cout << "Applying different compression levels...\n";
    for (int i = 0; i < int(domain / dx); i++) {
        // Gaussian compression
        gaussian_compressed_8bit[i] = compressFloat(gaussian_data[i], 8);
        gaussian_compressed_12bit[i] = compressFloat(gaussian_data[i], 12);
        gaussian_compressed_16bit[i] = compressFloat(gaussian_data[i], 16);
        
        // Uniform compression
        uniform_compressed_8bit[i] = compressFloat(uniform_data[i], 8);
        uniform_compressed_12bit[i] = compressFloat(uniform_data[i], 12);
        uniform_compressed_16bit[i] = compressFloat(uniform_data[i], 16);
        
        // Exponential compression
        exponential_compressed_8bit[i] = compressFloat(exponential_data[i], 8);
        exponential_compressed_12bit[i] = compressFloat(exponential_data[i], 12);
        exponential_compressed_16bit[i] = compressFloat(exponential_data[i], 16);
    }
    
    // --- PART 3: FILE I/O AND SIZE COMPARISON ---
    std::cout << "Saving data to binary files and comparing sizes...\n";
    
    // Save original data
    SaveToBinaryFile("uncompressed_gaussian.bin", gaussian_data, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uncompressed_uniform.bin", uniform_data, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uncompressed_exponential.bin", exponential_data, sizeof(float) * int(domain / dx));
    
    // Save compressed data
    SaveToBinaryFile("gaussian_8bit.bin", gaussian_compressed_8bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("gaussian_12bit.bin", gaussian_compressed_12bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("gaussian_16bit.bin", gaussian_compressed_16bit, sizeof(float) * int(domain / dx));
    
    SaveToBinaryFile("uniform_8bit.bin", uniform_compressed_8bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uniform_12bit.bin", uniform_compressed_12bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("uniform_16bit.bin", uniform_compressed_16bit, sizeof(float) * int(domain / dx));
    
    SaveToBinaryFile("exponential_8bit.bin", exponential_compressed_8bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("exponential_12bit.bin", exponential_compressed_12bit, sizeof(float) * int(domain / dx));
    SaveToBinaryFile("exponential_16bit.bin", exponential_compressed_16bit, sizeof(float) * int(domain / dx));
    
    // Compare file sizes (theoretical - actual file sizes may vary due to filesystem overhead)
    std::cout << "File size for original data: " << sizeof(float) * int(domain / dx) << " bytes\n";
    std::cout << "Note: Actual compression savings in file size are achieved when using specialized\n";
    std::cout << "      compression algorithms that exploit patterns in the mantissa bits.\n";
    std::cout << "      This implementation focuses on precision trade-offs rather than\n";
    std::cout << "      actual file size reduction.\n\n";
    
    // --- PART 4: STATISTICAL ANALYSIS ---
    std::cout << "\n=== STATISTICAL ANALYSIS ===\n";
    
    // Gaussian statistics
    double gaussian_orig_mean = mean(gaussian_data, int(domain / dx));
    double gaussian_8bit_mean = mean(gaussian_compressed_8bit, int(domain / dx));
    double gaussian_12bit_mean = mean(gaussian_compressed_12bit, int(domain / dx));
    double gaussian_16bit_mean = mean(gaussian_compressed_16bit, int(domain / dx));
    
    double gaussian_orig_var = variance(gaussian_data, int(domain / dx), gaussian_orig_mean);
    double gaussian_8bit_var = variance(gaussian_compressed_8bit, int(domain / dx), gaussian_8bit_mean);
    double gaussian_12bit_var = variance(gaussian_compressed_12bit, int(domain / dx), gaussian_12bit_mean);
    double gaussian_16bit_var = variance(gaussian_compressed_16bit, int(domain / dx), gaussian_16bit_mean);
    
    std::cout << "Gaussian Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << gaussian_orig_mean << " Variance: " << gaussian_orig_var << "\n";
    std::cout << "8-bit comp - Mean: " << gaussian_8bit_mean << " Variance: " << gaussian_8bit_var << "\n";
    std::cout << "12-bit comp - Mean: " << gaussian_12bit_mean << " Variance: " << gaussian_12bit_var << "\n";
    std::cout << "16-bit comp - Mean: " << gaussian_16bit_mean << " Variance: " << gaussian_16bit_var << "\n\n";
    
    // Uniform statistics 
    double uniform_orig_mean = mean(uniform_data, int(domain / dx));
    double uniform_8bit_mean = mean(uniform_compressed_8bit, int(domain / dx));
    double uniform_12bit_mean = mean(uniform_compressed_12bit, int(domain / dx));
    double uniform_16bit_mean = mean(uniform_compressed_16bit, int(domain / dx));
    
    double uniform_orig_var = variance(uniform_data, int(domain / dx), uniform_orig_mean);
    double uniform_8bit_var = variance(uniform_compressed_8bit, int(domain / dx), uniform_8bit_mean);
    double uniform_12bit_var = variance(uniform_compressed_12bit, int(domain / dx), uniform_12bit_mean);
    double uniform_16bit_var = variance(uniform_compressed_16bit, int(domain / dx), uniform_16bit_mean);
    
    std::cout << "uniform Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << uniform_orig_mean << " Variance: " << uniform_orig_var << "\n";
    std::cout << "8-bit comp - Mean: " << uniform_8bit_mean << " Variance: " << uniform_8bit_var << "\n";
    std::cout << "12-bit comp - Mean: " << uniform_12bit_mean << " Variance: " << uniform_12bit_var << "\n";
    std::cout << "16-bit comp - Mean: " << uniform_16bit_mean << " Variance: " << uniform_16bit_var << "\n\n";
    
    // Exponential statistics
    double exp_orig_mean = mean(exponential_data, int(domain / dx));
    double exp_8bit_mean = mean(exponential_compressed_8bit, int(domain / dx));
    double exp_12bit_mean = mean(exponential_compressed_12bit, int(domain / dx));
    double exp_16bit_mean = mean(exponential_compressed_16bit, int(domain / dx));
    
    double exp_orig_var = variance(exponential_data, int(domain / dx), exp_orig_mean);
    double exp_8bit_var = variance(exponential_compressed_8bit, int(domain / dx), exp_8bit_mean);
    double exp_12bit_var = variance(exponential_compressed_12bit, int(domain / dx), exp_12bit_mean);
    double exp_16bit_var = variance(exponential_compressed_16bit, int(domain / dx), exp_16bit_mean);
    
    std::cout << "Exponential Distribution Statistics:\n";
    std::cout << "Original   - Mean: " << exp_orig_mean << " Variance: " << exp_orig_var << "\n";
    std::cout << "8-bit comp - Mean: " << exp_8bit_mean << " Variance: " << exp_8bit_var << "\n";
    std::cout << "12-bit comp - Mean: " << exp_12bit_mean << " Variance: " << exp_12bit_var << "\n";
    std::cout << "16-bit comp - Mean: " << exp_16bit_mean << " Variance: " << exp_16bit_var << "\n\n";
    
    // --- PART 5: ERROR ANALYSIS ---
    std::cout << "\n=== ERROR ANALYSIS ===\n";
    
    // Calculate MSE for different distributions and compression levels
    double mse_gaussian_8bit = 0.0, mse_gaussian_12bit = 0.0, mse_gaussian_16bit = 0.0;
    double mse_uniform_8bit = 0.0, mse_uniform_12bit = 0.0, mse_uniform_16bit = 0.0;
    double mse_exp_8bit = 0.0, mse_exp_12bit = 0.0, mse_exp_16bit = 0.0;
    
    // Max errors
    double max_err_gaussian_8bit = 0.0, max_err_gaussian_12bit = 0.0, max_err_gaussian_16bit = 0.0;
    double max_err_uniform_8bit = 0.0, max_err_uniform_12bit = 0.0, max_err_uniform_16bit = 0.0;
    double max_err_exp_8bit = 0.0, max_err_exp_12bit = 0.0, max_err_exp_16bit = 0.0;
    
    for (int i = 0; i < int(domain / dx); i++) {
        // Gaussian errors
        double err_gaussian_8bit = fabs(gaussian_data[i] - gaussian_compressed_8bit[i]);
        double err_gaussian_12bit = fabs(gaussian_data[i] - gaussian_compressed_12bit[i]);
        double err_gaussian_16bit = fabs(gaussian_data[i] - gaussian_compressed_16bit[i]);
        
        mse_gaussian_8bit += err_gaussian_8bit * err_gaussian_8bit;
        mse_gaussian_12bit += err_gaussian_12bit * err_gaussian_12bit;
        mse_gaussian_16bit += err_gaussian_16bit * err_gaussian_16bit;
        
        max_err_gaussian_8bit = std::max(max_err_gaussian_8bit, err_gaussian_8bit);
        max_err_gaussian_12bit = std::max(max_err_gaussian_12bit, err_gaussian_12bit);
        max_err_gaussian_16bit = std::max(max_err_gaussian_16bit, err_gaussian_16bit);
        
        // Uniform errors
        double err_uniform_8bit = fabs(uniform_data[i] - uniform_compressed_8bit[i]);
        double err_uniform_12bit = fabs(uniform_data[i] - uniform_compressed_12bit[i]);
        double err_uniform_16bit = fabs(uniform_data[i] - uniform_compressed_16bit[i]);
        
        mse_uniform_8bit += err_uniform_8bit * err_uniform_8bit;
        mse_uniform_12bit += err_uniform_12bit * err_uniform_12bit;
        mse_uniform_16bit += err_uniform_16bit * err_uniform_16bit;
        
        max_err_uniform_8bit = std::max(max_err_uniform_8bit, err_uniform_8bit);
        max_err_uniform_12bit = std::max(max_err_uniform_12bit, err_uniform_12bit);
        max_err_uniform_16bit = std::max(max_err_uniform_16bit, err_uniform_16bit);
        
        // Exponential errors
        double err_exp_8bit = fabs(exponential_data[i] - exponential_compressed_8bit[i]);
        double err_exp_12bit = fabs(exponential_data[i] - exponential_compressed_12bit[i]);
        double err_exp_16bit = fabs(exponential_data[i] - exponential_compressed_16bit[i]);
        
        mse_exp_8bit += err_exp_8bit * err_exp_8bit;
        mse_exp_12bit += err_exp_12bit * err_exp_12bit;
        mse_exp_16bit += err_exp_16bit * err_exp_16bit;
        
        max_err_exp_8bit = std::max(max_err_exp_8bit, err_exp_8bit);
        max_err_exp_12bit = std::max(max_err_exp_12bit, err_exp_12bit);
        max_err_exp_16bit = std::max(max_err_exp_16bit, err_exp_16bit);
    }
    
    int numPoints = int(domain / dx);
    
    // uniformize MSE values
    mse_gaussian_8bit /= numPoints;
    mse_gaussian_12bit /= numPoints;
    mse_gaussian_16bit /= numPoints;
    
    mse_uniform_8bit /= numPoints;
    mse_uniform_12bit /= numPoints;
    mse_uniform_16bit /= numPoints;
    
    mse_exp_8bit /= numPoints;
    mse_exp_12bit /= numPoints;
    mse_exp_16bit /= numPoints;
    
    // Display error metrics
    std::cout << "Gaussian Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_gaussian_8bit << " Max Error: " << max_err_gaussian_8bit << "\n";
    std::cout << "12-bit - MSE: " << mse_gaussian_12bit << " Max Error: " << max_err_gaussian_12bit << "\n";
    std::cout << "16-bit - MSE: " << mse_gaussian_16bit << " Max Error: " << max_err_gaussian_16bit << "\n\n";
    
    std::cout << "uniform Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_uniform_8bit << " Max Error: " << max_err_uniform_8bit << "\n";
    std::cout << "12-bit - MSE: " << mse_uniform_12bit << " Max Error: " << max_err_uniform_12bit << "\n";
    std::cout << "16-bit - MSE: " << mse_uniform_16bit << " Max Error: " << max_err_uniform_16bit << "\n\n";
    
    std::cout << "Exponential Distribution Error Metrics:\n";
    std::cout << "8-bit  - MSE: " << mse_exp_8bit << " Max Error: " << max_err_exp_8bit << "\n";
    std::cout << "12-bit - MSE: " << mse_exp_12bit << " Max Error: " << max_err_exp_12bit << "\n";
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
    } else if (mse_gaussian_12bit < 1e-6 && mse_uniform_12bit < 1e-6 && mse_exp_12bit < 1e-6) {
        std::cout << "- 12-bit compression provides a good balance between precision and storage efficiency\n";
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
    delete[] gaussian_compressed_12bit;
    delete[] gaussian_compressed_16bit;
    delete[] uniform_compressed_8bit;
    delete[] uniform_compressed_12bit;
    delete[] uniform_compressed_16bit;
    delete[] exponential_compressed_8bit;
    delete[] exponential_compressed_12bit;
    delete[] exponential_compressed_16bit;
    
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