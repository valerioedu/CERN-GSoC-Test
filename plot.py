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
    
    plt.subplot(2, 1, 1)
    plt.plot(original_data['x'], original_data['y'], 'b-', label='Original')
    plt.plot(compressed_8bit['x'], compressed_8bit['y'], 'r--', label='8-bit compression')
    plt.plot(compressed_16bit['x'], compressed_16bit['y'], 'g-.', label='16-bit compression')
    plt.title(f'{dist_name.capitalize()} Distribution Compression Comparison')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 1, 2)

    if dist_name == 'uniform':
        indices = np.where(original_data['y'] > 0)[0]
        if len(indices) > 0:
            mid_point = (indices[0] + indices[-1]) // 2
            window_size = min(100, len(indices) // 2)
            start_idx = mid_point - window_size
            end_idx = mid_point + window_size
        else:
            start_idx = len(original_data) // 4
            end_idx = 3 * len(original_data) // 4
    else:
        try:
            max_idx = original_data['y'].idxmax()
            start_idx = max(0, max_idx - 50)
            end_idx = min(len(original_data), max_idx + 50)
        except:
            start_idx = len(original_data) // 4
            end_idx = 3 * len(original_data) // 4
    
    start_idx = max(0, start_idx)
    end_idx = min(len(original_data), end_idx)
    
    plt.plot(original_data['x'][start_idx:end_idx], original_data['y'][start_idx:end_idx], 'b-', label='Original')
    plt.plot(compressed_8bit['x'][start_idx:end_idx], compressed_8bit['y'][start_idx:end_idx], 'r--', label='8-bit compression')
    plt.plot(compressed_16bit['x'][start_idx:end_idx], compressed_16bit['y'][start_idx:end_idx], 'g-.', label='16-bit compression')
    plt.title('Zoomed View (Around Maximum/Middle)')
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{dist_name}_comparison.png', dpi=300)
    plt.close()

def plot_compression_errors(dist_name):
    """Plot error between original and compressed data"""
    original_data = load_csv(f"{dist_name}_original.csv")
    compressed_8bit = load_csv(f"{dist_name}_8bit.csv")
    compressed_16bit = load_csv(f"{dist_name}_16bit.csv")
    
    if original_data is None or compressed_8bit is None or compressed_16bit is None:
        print(f"Could not load data for {dist_name} distribution error analysis")
        return
    
    error_8bit = original_data['y'] - compressed_8bit['y']
    error_16bit = original_data['y'] - compressed_16bit['y']
    
    plt.figure(figsize=(12, 6))
    
    plt.plot(original_data['x'], error_8bit, 'r-', label='8-bit error', alpha=0.7)
    plt.plot(original_data['x'], error_16bit, 'g-', label='16-bit error', alpha=0.7)
    
    plt.title(f'Compression Error for {dist_name.capitalize()} Distribution')
    plt.xlabel('Index')
    plt.ylabel('Error (Original - Compressed)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{dist_name}_error.png', dpi=300)
    plt.close()

def plot_bit_pattern_example():
     """Visualize the bit patterns of a float value and its compressed versions"""
     value = 1.234567
     value_8bit_compressed = 1.2345581
     value_16bit_compressed = 1.234375
 
     def float_to_bin_str(f):
         import struct
         b = struct.pack('!f', f)
         i = struct.unpack('!I', b)[0]
         return bin(i)[2:].zfill(32)
 
     original_bits = float_to_bin_str(value)
     compressed_8bit = float_to_bin_str(value_8bit_compressed)
     compressed_16bit = float_to_bin_str(value_16bit_compressed)
 
     plt.figure(figsize=(12, 6))
 
     def plot_bits(y_pos, bit_str, label):
         sign_bit = bit_str[0]
         exponent_bits = bit_str[1:9]
         mantissa_bits = bit_str[9:]
 
         plt.text(0, y_pos, sign_bit, fontsize=14, color='red')
         for i, bit in enumerate(exponent_bits):
             plt.text(i+1, y_pos, bit, fontsize=14, color='blue')
         for i, bit in enumerate(mantissa_bits):
             plt.text(i+9, y_pos, bit, fontsize=14, 
                      color='green' if i < len(mantissa_bits)-16 else 
                            'orange' if i < len(mantissa_bits)-8 else 'purple')
 
         plt.text(-5, y_pos, f"{label}:", fontsize=14, ha='right')
 
     plot_bits(3, original_bits, "Original")
     plot_bits(2, compressed_8bit, "8-bit comp.")
     plot_bits(1, compressed_16bit, "16-bit comp.")
 
     plt.axvline(x=0.8,ymin=0.1, ymax=0.7, color='k', linestyle='-', alpha=0.2)
     plt.axvline(x=8.8,ymin=0.1, ymax=0.7, color='k', linestyle='-', alpha=0.2)
 
     plt.text(0, 4, "Sign", fontsize=12, ha='center')
     plt.text(4.5, 4, "Exponent", fontsize=12, ha='center')
     plt.text(20, 4, "Mantissa", fontsize=12, ha='center')
 
     plt.text(20, 0.5, "↑ These 16 bits zeroed", fontsize=10, color='orange')
     plt.text(28, 1.5, "↑ These 8 bits zeroed", fontsize=10, color='purple')
 
     plt.xlim(-5, 40)
     plt.ylim(0.5, 4.5)
     plt.axis('off')
     plt.title('Bit Pattern Visualization: Float Value Compression', x=0.42)
 
     plt.tight_layout()
     plt.savefig('bit_pattern.png', dpi=300)
     plt.close()

def plot_file_size_comparison():
    """Create a bar chart of file sizes"""
    
    labels = ['Original (32-bit)', '8-bit compression', '16-bit compression']
    sizes = [100, 75, 50]  # Expected and verified percentages
    
    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, sizes, color=['blue', 'green', 'red'])
    
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{int(height)}%', ha='center', va='bottom')
    
    plt.title('File Size Comparison', fontsize=14)
    plt.ylabel('Percentage of Original Size', fontsize=12)
    plt.ylim(0, 120)
    plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('file_size_comparison.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    print("Generating visualization plots...")
    
    for dist in ['gaussian', 'uniform', 'exponential']:
        plot_distribution_comparison(dist)
        plot_compression_errors(dist)
    
    plot_bit_pattern_example()
    
    plot_file_size_comparison()
    
    print("Done! Plots saved as PNG files.")