import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_data(file_path):
    # Read the data from the file into a DataFrame
    data = pd.read_csv(file_path, header=None, names=['event', 'time', 'value'])

    # Separate data for new_current, new_incumbent, and tenure events
    current_data = data[data['event'] == 'new_current']
    incumbent_data = data[data['event'] == 'new_incumbent']
    tenure_data = data[data['event'] == 'tenure']

    # Plot the data for new_current and new_incumbent events
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(current_data['time'], current_data['value'], marker='o', linestyle='-', label='current', markersize=0)
    ax1.plot(incumbent_data['time'], incumbent_data['value'], marker='x', linestyle='--', label='incumbent')

    # Add labels and title for the first y-axis
    ax1.set_xlabel('iteration')
    # ax1.set_ylabel('Value (new_current, new_incumbent)')
    # ax1.set_title('Plot of Data')

    # Plot the data for the tenure event on a secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(tenure_data['time'], tenure_data['value'], marker='s', linestyle='-.', color='m', label='tenure', markersize=0)
    # ax2.set_ylabel('Value (tenure)')

    # Add legend for both y-axes
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    # Show plot
    plt.grid(True)

    plt.savefig("out.pdf")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    plot_data(file_path)
