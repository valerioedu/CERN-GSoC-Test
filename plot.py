import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def load_csv(filename):
    """Load data from CSV file with x,y format"""
    if not os.path.exists(filename):
        print(f"Warning: File {filename} not found")
        return None
    
    data = pd.read_csv(filename, header=None, names=['x', 'y'])
    return data

def plot_distribution_comparison(dist_name):
    """Plot original vs compressed data for a specific distribution"""
    original_data = load_csv(f"{dist_name}_original.csv")
    compressed_8bit = load_csv(f"{dist_name}_8bit.csv")
    compressed_16bit = load_csv(f"{dist_name}_16bit.csv")
    
    if original_data is None or compressed_8bit is None or compressed_16bit is None:
        print(f"Could not load data for {dist_name} distribution")
        return
    
    plt.figure(figsize=(12, 8))
    
    # Main plot - full range
    plt.subplot(2, 1, 1)
    plt.plot(original_data['x'], original_data['y'], 'b-', label='Original')
    plt.plot(compressed_8bit['x'], compressed_8bit['y'], 'r--', label='8-bit compression')
    plt.plot(compressed_16bit['x'], compressed_16bit['y'], 'g-.', label='16-bit compression')
    plt.title(f'{dist_name.capitalize()} Distribution Compression Comparison')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Zoomed plot - to see differences better
    plt.subplot(2, 1, 2)
    
    # Find index of maximum value to focus plot around it
    max_idx = original_data['y'].idxmax()
    start_idx = max(0, max_idx - 50)
    end_idx = min(len(original_data), max_idx + 50)
    
    plt.plot(original_data['x'][start_idx:end_idx], original_data['y'][start_idx:end_idx], 'b-', label='Original')
    plt.plot(compressed_8bit['x'][start_idx:end_idx], compressed_8bit['y'][start_idx:end_idx], 'r--', label='8-bit compression')
    plt.plot(compressed_16bit['x'][start_idx:end_idx], compressed_16bit['y'][start_idx:end_idx], 'g-.', label='16-bit compression')
    plt.title('Zoomed View (Around Maximum)')
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{dist_name}_comparison.png', dpi=300)
    plt.close()

def plot_error_comparison():
    """Plot error metrics across distributions and compression levels"""
    distributions = ['gaussian', 'uniform', 'exponential']
    compression_levels = ['8bit', '12bit', '16bit']
    
    # Sample data - replace with actual error metrics from your output
    # You could extract these from the binary files or hardcode from the console output
    mse_data = {
        'gaussian': [1e-7, 1e-6, 1e-5],  # MSE for 8bit, 12bit, 16bit
        'uniform': [1.5e-7, 1.5e-6, 1.5e-5],
        'exponential': [2e-7, 2e-6, 2e-5]
    }
    
    max_error_data = {
        'gaussian': [0.0001, 0.001, 0.01],  # Max errors for 8bit, 12bit, 16bit
        'uniform': [0.00015, 0.0015, 0.015],
        'exponential': [0.0002, 0.002, 0.02]
    }
    
    # Prepare bar chart data
    plt.figure(figsize=(14, 10))
    
    # MSE comparison
    plt.subplot(2, 1, 1)
    x = np.arange(len(distributions))
    width = 0.25
    
    plt.bar(x - width, [mse_data[d][0] for d in distributions], width, label='8-bit')
    plt.bar(x, [mse_data[d][1] for d in distributions], width, label='12-bit')
    plt.bar(x + width, [mse_data[d][2] for d in distributions], width, label='16-bit')
    
    plt.yscale('log')  # Use log scale for better visibility of differences
    plt.title('Mean Squared Error by Distribution and Compression Level')
    plt.xticks(x, distributions)
    plt.ylabel('Mean Squared Error (log scale)')
    plt.legend()
    plt.grid(True, axis='y', alpha=0.3)
    
    # Max Error comparison
    plt.subplot(2, 1, 2)
    
    plt.bar(x - width, [max_error_data[d][0] for d in distributions], width, label='8-bit')
    plt.bar(x, [max_error_data[d][1] for d in distributions], width, label='12-bit')
    plt.bar(x + width, [max_error_data[d][2] for d in distributions], width, label='16-bit')
    
    plt.yscale('log')  # Use log scale for better visibility of differences
    plt.title('Maximum Error by Distribution and Compression Level')
    plt.xticks(x, distributions)
    plt.ylabel('Maximum Error (log scale)')
    plt.legend()
    plt.grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('error_comparison.png', dpi=300)
    plt.close()

def plot_bit_pattern_example():
    """Visualize the bit patterns of a float value and its compressed versions"""
    # Original value and compressed versions for demonstration
    value = 1.234567
    value_8bit_compressed = 1.2345581
    value_16bit_compressed = 1.234375
    
    # Convert to binary representation (just for demonstration)
    def float_to_bin_str(f):
        # This is simplified - in reality you'd need proper IEEE-754 conversion
        import struct
        b = struct.pack('!f', f)
        i = struct.unpack('!I', b)[0]
        return bin(i)[2:].zfill(32)
    
    original_bits = float_to_bin_str(value)
    compressed_8bit = float_to_bin_str(value_8bit_compressed)
    compressed_16bit = float_to_bin_str(value_16bit_compressed)
    
    # Create a visual representation
    plt.figure(figsize=(12, 6))
    
    # Helper to display bits with coloring
    def plot_bits(y_pos, bit_str, label):
        sign_bit = bit_str[0]
        exponent_bits = bit_str[1:9]
        mantissa_bits = bit_str[9:]
        
        # Plot each section with different colors
        plt.text(0, y_pos, sign_bit, fontsize=14, color='red')
        for i, bit in enumerate(exponent_bits):
            plt.text(i+1, y_pos, bit, fontsize=14, color='blue')
        for i, bit in enumerate(mantissa_bits):
            plt.text(i+9, y_pos, bit, fontsize=14, 
                     color='green' if i < len(mantissa_bits)-16 else 
                           'orange' if i < len(mantissa_bits)-8 else 'purple')
        
        # Add label
        plt.text(-5, y_pos, f"{label}:", fontsize=14, ha='right')
    
    plot_bits(3, original_bits, "Original")
    plot_bits(2, compressed_8bit, "8-bit comp.")
    plot_bits(1, compressed_16bit, "16-bit comp.")
    
    # Add vertical separators and labels
    plt.axvline(x=0.8,ymin=0.1, ymax=0.7, color='k', linestyle='-', alpha=0.2)
    plt.axvline(x=8.8,ymin=0.1, ymax=0.7, color='k', linestyle='-', alpha=0.2)
    
    plt.text(0, 4, "Sign", fontsize=12, ha='center')
    plt.text(4.5, 4, "Exponent", fontsize=12, ha='center')
    plt.text(20, 4, "Mantissa", fontsize=12, ha='center')
    
    # Mark the bits that are zeroed in each compression
    plt.text(20, 0.5, "↑ These 16 bits zeroed", fontsize=10, color='orange')
    plt.text(28, 1.5, "↑ These 8 bits zeroed", fontsize=10, color='purple')
    
    plt.xlim(-5, 40)
    plt.ylim(0.5, 4.5)
    plt.axis('off')
    plt.title('Bit Pattern Visualization: Float Value Compression')
    
    plt.tight_layout()
    plt.savefig('bit_pattern.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    print("Generating visualization plots...")
    
    # Plot distribution comparisons
    for dist in ['gaussian', 'uniform', 'exponential']:
        plot_distribution_comparison(dist)
    
    # Plot error metrics comparison
    plot_error_comparison()
    
    # Plot bit pattern example
    plot_bit_pattern_example()
    
    print("Done! Plots saved as PNG files.")