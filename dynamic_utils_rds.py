"""
Enhanced utility functions for column-based dynamic options with column names
"""
import sys
import os
import subprocess

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

        # Find Rscript in Galaxy's environment
        import shutil
        rscript_path = None
        
        # Try multiple approaches to find Rscript
        rscript_candidates = [
            'Rscript',  # Try PATH first
            '/usr/local/bin/Rscript',
            '/usr/bin/Rscript',
            '/opt/conda/bin/Rscript',  # Common in Galaxy conda environments
            '/galaxy/server/database/dependencies/_conda/bin/Rscript',  # Galaxy conda
            '/tool_deps/_conda/bin/Rscript',  # Tool dependencies
            '/conda-envs/*/bin/Rscript',  # Galaxy conda environments pattern
        ]
        
        # Also check environment variables that Galaxy might set
        if 'CONDA_PREFIX' in os.environ:
            rscript_candidates.insert(1, os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'Rscript'))
        if 'R_HOME' in os.environ:
            rscript_candidates.insert(1, os.path.join(os.environ['R_HOME'], 'bin', 'Rscript'))
        
        # Check for conda environments in common Galaxy locations
        import glob
        for conda_env_pattern in ['/opt/conda/envs/*/bin/Rscript', '/conda-envs/*/bin/Rscript']:
            rscript_candidates.extend(glob.glob(conda_env_pattern))
        
        for candidate in rscript_candidates:
            if shutil.which(candidate) or (os.path.isabs(candidate) and os.path.isfile(candidate)):
                rscript_path = candidate
                break
        
        if not rscript_path:
            # As a last resort, try to find R and construct Rscript path
            r_path = shutil.which('R')
            if r_path:
                potential_rscript = os.path.join(os.path.dirname(r_path), 'Rscript')
                if os.path.isfile(potential_rscript):
                    rscript_path = potential_rscript
        
        if not rscript_path:
            env_info = f"PATH: {os.environ.get('PATH', 'Not set')}\nCONDA_PREFIX: {os.environ.get('CONDA_PREFIX', 'Not set')}\nR_HOME: {os.environ.get('R_HOME', 'Not set')}"
            return [(f"Error: Rscript not found in Galaxy environment. Environment: {env_info}", "error", True)]
        
        # Get the directory where this script is located to find Get_metadata.R
        # Try multiple approaches to find the script directory
        script_dir = None
        
        # Method 1: Use __tool_directory__ if available (Galaxy sets this)
        try:
            if '__tool_directory__' in globals():
                script_dir = globals()['__tool_directory__']
        except:
            pass
        
        # Method 2: Use __file__ if available
        if not script_dir:
            try:
                if '__file__' in globals():
                    script_dir = os.path.dirname(os.path.abspath(globals()['__file__']))
            except:
                pass
        
        # Method 3: Try to find Get_metadata.R in common locations
        if not script_dir:
            potential_dirs = [
                '.',  # Current working directory
                os.getcwd(),  # Explicitly get current working directory
                '/galaxy_tools',  # Common Galaxy tool directory
                '/tool_deps',  # Tool dependencies directory
                os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else None,
            ]
            
            # Remove None values
            potential_dirs = [d for d in potential_dirs if d is not None]
            
            # Also try relative to where dynamic_utils_rds.py might be
            try:
                import dynamic_utils_rds
                if hasattr(dynamic_utils_rds, '__file__'):
                    potential_dirs.append(os.path.dirname(os.path.abspath(dynamic_utils_rds.__file__)))
            except:
                pass
            
            for potential_dir in potential_dirs:
                try:
                    if os.path.exists(os.path.join(potential_dir, "Get_metadata.R")):
                        script_dir = potential_dir
                        break
                except:
                    continue
        
        if not script_dir:
            script_dir = '.'  # Fallback to current directory
        
        # Try both with and without .R extension in case of confusion
        possible_r_script_names = ["Get_metadata.R", "Get_metadata"]
        r_script_path = None
        
        for script_name in possible_r_script_names:
            test_path = os.path.join(script_dir, script_name)
            if os.path.exists(test_path):
                r_script_path = test_path
                break
        
        # Enhanced error message if script not found
        if not r_script_path:
            # List files in the script directory for debugging
            try:
                files_in_dir = os.listdir(script_dir)
                r_files = [f for f in files_in_dir if f.endswith('.R')]
                files_info = f"Working dir: {os.getcwd()}, Script dir: {script_dir}, R files found: {r_files}, All files: {files_in_dir[:10]}"
            except Exception as e:
                files_info = f"Could not list files in {script_dir}: {e}"
            
            return [(f"Error: Get_metadata.R not found. {files_info}", "error", True)]
        
        # Test if Seurat is available in this R environment
        test_cmd = f'"{rscript_path}" -e "library(Seurat); cat(\\"Seurat OK\\")"'
        try:
            test_result = subprocess.run(test_cmd, shell=True, check=True, capture_output=True, text=True)
            if "Seurat OK" not in test_result.stdout:
                return [(f"Error: Seurat not available in R environment. Output: {test_result.stdout}", "error", True)]
        except subprocess.CalledProcessError as e:
            return [(f"Error: Seurat library test failed. Stderr: {e.stderr}", "error", True)]

        r_cmd = f'"{rscript_path}" --vanilla "{r_script_path}" "{file_path}"'

        # Debug information - log the command and environment
        debug_msg = f"""
=== R Script Execution Debug ===
Command: {r_cmd}
Working directory: {os.getcwd()}
Script directory detected: {script_dir}
Script path: {r_script_path}
Script exists: {os.path.exists(r_script_path)}
Input file exists: {os.path.exists(file_path)}
Available globals: {list(globals().keys())}
Environment:
  PATH: {os.environ.get('PATH', 'Not set')}
  CONDA_PREFIX: {os.environ.get('CONDA_PREFIX', 'Not set')}
  R_HOME: {os.environ.get('R_HOME', 'Not set')}
"""
        
        try:
            with open('galaxy_dynamic_debug.log', 'a') as f:
                f.write(debug_msg)
        except:
            pass

        try:
            result = subprocess.run(r_cmd, shell=True, check=True, capture_output=True, text=True, 
                                  cwd=script_dir)  # Run from script directory
        except subprocess.CalledProcessError as e:
            error_details = f"""R Script Error:
Command: {r_cmd}
Return code: {e.returncode}
Stdout: {e.stdout}
Stderr: {e.stderr}
Working dir: {script_dir}
"""
            return [(f"Error running R script: {error_details}", "error", True)]
        
        # Check if the metadata file was created
        if not os.path.exists("combined_metadata.csv"):
            return [("Error: metadata file not created", "error", True)]
        
        with open("combined_metadata.csv", 'r', encoding='utf-8', errors='replace') as f:
            first_line = f.readline().strip()
            
            if not first_line:
                return [("Empty file", "empty", True)]
            
            # Check if the file contains an error message
            if first_line.startswith("Error"):
                return [(first_line, "error", True)]
            
            # Split by comma (since we know it's a CSV file)
            headers = first_line.split(',')
            
            # Create options (column_name, column_name, selected)
            # Store just the column name as the value, we'll look up the index later
            for i, header in enumerate(headers):
                header = header.strip().strip('"')  # Remove quotes and whitespace
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
        # Check if metadata file exists (should have been created by get_column_names_options)
        if not os.path.exists("combined_metadata.csv"):
            return [("Metadata file not found. Please select a column first.", "error", True)]
        
        # First, get the column headers to find the index of the selected column
        file_path = dataset.get_file_name() if hasattr(dataset, 'get_file_name') else dataset.file_name
        
        with open("combined_metadata.csv", 'r', encoding='utf-8', errors='replace') as f:
            first_line = f.readline().strip()
            
            if not first_line:
                return [("Empty file", "empty", True)]
            
            # Split by comma (since we know it's a CSV file)
            headers = first_line.split(',')
            headers = [h.strip().strip('"') for h in headers]  # Remove quotes and whitespace
        
        # Find the column index based on the column name
        col_index = -1
        for i, header in enumerate(headers):
            if header == selected_column:
                col_index = i
                break
        
        if col_index == -1:
            return [("Column not found", "error", True)]
        
        # Now read the data and get unique values from that column
        # Read from the CSV metadata file, not the original RDS file
        unique_values = set()
        
        with open("combined_metadata.csv", 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
            
            # Skip header (first line)
            data_lines = lines[1:] if len(lines) > 1 else []
            
            for line in data_lines:
                line = line.strip()
                if line:
                    # Split by comma (since we know it's a CSV file)
                    columns = line.split(',')
                    if col_index < len(columns):
                        value = columns[col_index].strip()
                        # Remove quotes if present
                        if value.startswith('"') and value.endswith('"'):
                            value = value[1:-1]
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
