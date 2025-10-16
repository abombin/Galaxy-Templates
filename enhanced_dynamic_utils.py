"""
Enhanced utility functions for column-based dynamic options with column names
"""
import sys
import os

# Print Python version and executable path for debugging
debug_info = f"""Python version: {sys.version}
Python executable: {sys.executable}
Current working directory: {os.getcwd()}
"""

print(debug_info)
print(debug_info, file=sys.stderr)

# Also write to a debug file if possible
try:
    with open('galaxy_dynamic_debug.log', 'w') as f:
        f.write(f"=== Dynamic Options Debug - {os.getpid()} ===\n")
        f.write(debug_info)
except:
    pass  # Ignore if we can't write to the file

def get_column_names_options(dataset):
    """
    Get column names from the dataset header as options.
    """
    options = []
    
    if dataset is None:
        return [("No dataset selected", "none", True)]
    
    try:
        # Read the dataset file to get column headers
        file_path = dataset.get_file_name() if hasattr(dataset, 'get_file_name') else dataset.file_name
        
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            
            if not first_line:
                return [("Empty file", "empty", True)]
            
            # Split by tab or comma depending on file format
            headers = first_line.split('\t') if '\t' in first_line else first_line.split(',')
            
            # Create options (column_name, column_name, selected)
            # Store just the column name as the value, we'll look up the index later
            for i, header in enumerate(headers):
                header = header.strip()
                if header:
                    options.append((header, header, i == 0))
                else:
                    col_name = f"Column_{i+1}"
                    options.append((col_name, col_name, i == 0))
                    
    except Exception as e:
        return [("Error reading headers: " + str(e), "error", True)]
    
    return options

def get_column_values_options(dataset, selected_column=None):
    """
    Get unique values from the selected column as options.
    """
    options = []
    
    if dataset is None or selected_column is None:
        return [("Select a column first", "none", True)]
    
    try:
        # First, get the column headers to find the index of the selected column
        file_path = dataset.get_file_name() if hasattr(dataset, 'get_file_name') else dataset.file_name
        
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            
            if not first_line:
                return [("Empty file", "empty", True)]
            
            # Split by tab or comma depending on file format
            headers = first_line.split('\t') if '\t' in first_line else first_line.split(',')
            headers = [h.strip() for h in headers]
        
        # Find the column index based on the column name
        col_index = -1
        for i, header in enumerate(headers):
            if header == selected_column:
                col_index = i
                break
        
        if col_index == -1:
            return [("Column not found", "error", True)]
        
        # Now read the data and get unique values from that column
        unique_values = set()
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
            # Skip header (first line)
            data_lines = lines[1:] if len(lines) > 1 else []
            
            for line in data_lines:
                line = line.strip()
                if line:
                    # Split by tab or comma depending on file format
                    columns = line.split('\t') if '\t' in line else line.split(',')
                    if col_index < len(columns):
                        value = columns[col_index].strip()
                        if value and value != '':  # Only add non-empty values
                            unique_values.add(value)
        
        # Convert to sorted list and create options
        sorted_values = sorted(list(unique_values))
        
        if not sorted_values:
            return [(f"No values found in '{selected_column}'", "empty", True)]
        
        # Create options (value, label, selected)
        # Limit to first 50 unique values to avoid overwhelming the UI
        limited_values = sorted_values[:50]
        for i, value in enumerate(limited_values):
            options.append((value, value, i == 0))
            
        # Add note if there are more values
        if len(sorted_values) > 50:
            options.append((f"... and {len(sorted_values) - 50} more values", "more", False))
            
    except Exception as e:
        return [("Error reading column: " + str(e), "error", True)]
    
    return options
