import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Optional ace_tools integration
try:
    import ace_tools as tools
    ace_tools_available = True
except ImportError:
    ace_tools_available = False

def process_csv_files(folder_path):
    """
    Processes CSV files in the specified folder to analyze CDR3 length distributions and top V-calls.

    Args:
        folder_path (str): Path to the folder containing CSV files.
    """
    analysis_folder = os.path.join(folder_path, "analysis")
    os.makedirs(analysis_folder, exist_ok=True)
    
    for file in os.listdir(folder_path):
        if file.endswith(".csv"):
            file_path = os.path.join(folder_path, file)
            df = pd.read_csv(file_path)
            
            # Extract filename components for naming
            base_name = file.split(".")[0]  # Everything before the first dot
            suffix = file.split("_")[-2]    # Last 2 chars before .csv
            output_name = f"{base_name}_{suffix}"
            
            # Extract length columns
            length_columns = [col for col in df.columns if col.isdigit()]
            
            # Normalize percentages to 0-100 scale and round to 1 decimal place
            df[length_columns] = (df[length_columns] * 100 / df[length_columns].sum().sum()).round(1)
            
            # Compute total percentage for each length and save as TSV
            length_distribution = df[length_columns].sum()
            length_tsv_path = os.path.join(analysis_folder, f"{output_name}_length_distribution.tsv")
            length_distribution.to_csv(length_tsv_path, sep='\t', header=['Percentage'])
            
            # Plot and save CDR3 length distribution as PNG
            _plot_length_distribution(length_distribution, output_name, analysis_folder)
            
            # Analyze top 7 V-calls by total percentage and save
            top_v_calls = _analyze_top_v_calls(df, length_columns)
            _save_top_v_calls(top_v_calls, output_name, analysis_folder)
            
            # Display table in ace_tools if available
            if ace_tools_available:
                tools.display_dataframe_to_user(
                    name=f"Top 7 V-Calls - {output_name}",
                    dataframe=top_v_calls[["v_call", "Total_Percentage", "Mean_Length"]]
                )
            
            print(f"Analysis complete for {file}! Files saved to: {analysis_folder}")


def _plot_length_distribution(length_distribution, output_name, analysis_folder):
    """Helper function to plot CDR3 length distribution."""
    plt.figure(figsize=(10, 6))
    plt.bar(length_distribution.index.astype(int), length_distribution.values)
    plt.xlabel("CDR3 Length")
    plt.ylabel("Percentage")
    plt.title(f"CDR3 Length Distribution - {output_name}")
    plt.grid(axis="y")
    
    graph_path = os.path.join(analysis_folder, f"{output_name}_cdr3_length_distribution.png")
    plt.savefig(graph_path)
    plt.close()


def _analyze_top_v_calls(df, length_columns):
    """Helper function to analyze top V-calls and compute mean lengths."""
    df["Total_Percentage"] = df[length_columns].sum(axis=1).round(2)
    top_v_calls = df.sort_values("Total_Percentage", ascending=False).head(7)
    
    mean_lengths = []
    for _, row in top_v_calls.iterrows():
        values = []
        for col in length_columns:
            count = int(row.get(col, 0))
            values.extend([int(col)] * count)
        mean_length = np.round(np.mean(values), 2) if values else 0
        mean_lengths.append(mean_length)
    
    top_v_calls["Mean_Length"] = mean_lengths
    return top_v_calls


def _save_top_v_calls(top_v_calls, output_name, analysis_folder):
    """Helper function to save the top V-calls data and visualization."""
    table_path = os.path.join(analysis_folder, f"{output_name}_top_7_v_calls.tsv")
    top_v_calls[["v_call", "Total_Percentage", "Mean_Length"]].to_csv(table_path, sep='\t', index=False)
    
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.axis('tight')
    ax.axis('off')
    
    column_labels = ["v_call", "Total_Percentage", "Mean_Length"]
    col_colors = ["#DDEBF7"] * len(column_labels)
    
    ax.table(
        cellText=top_v_calls[["v_call", "Total_Percentage", "Mean_Length"]].values,
        colLabels=column_labels,
        cellLoc='center',
        loc='center',
        colColours=col_colors
    )
    
    table_png_path = os.path.join(analysis_folder, f"{output_name}_top_7_v_calls.png")
    plt.savefig(table_png_path, bbox_inches='tight', dpi=300)
    plt.close()


# Allow running this script standalone for testing
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Error: No folder path provided. Usage: python3 vlength.py <folder_path>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    process_csv_files(folder_path)