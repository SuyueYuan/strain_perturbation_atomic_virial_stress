import os
import pandas as pd

# Constant for unit conversion from eV*A to GPa
evA2Gpa = 160.2

def read_lammps_data(filename, pe_column_index, vol_column_index):
    """
    Read a LAMMPS dump file, extracting atomic potential energy (pe) and atomic Voronoi volume (vol).
    Returns a DataFrame with columns ["pe", "vol"].
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ITEM: ATOMS"):
                for data_line in f:
                    columns = data_line.split()
                    data.append((
                        float(columns[pe_column_index]),
                        float(columns[vol_column_index])
                    ))
    return pd.DataFrame(data, columns=["pe", "vol"])


def compute_diff_virial(df_P, df_N, denom):
    """
    Compute difference in potential energy
    """
    return (df_P.pe - df_N.pe) * denom


def write_cfg_with_diff(original_file, output_file, diffs_dict):
    """
    Write out the dump.cfg file. 
    """
    directions_order = ["11","22","33","12","13","23"]
    diffs_lists = {k: diffs_dict[k].tolist() for k in directions_order}    
    n_atoms = len(diffs_lists["11"])
    
    with open(original_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines_buffer = []
        atom_index = 0

        for line_num, line in enumerate(f_in):
            if line_num < 8:
                # Copy header lines
                lines_buffer.append(line)
            elif line_num == 8:
                # Append new column headers
                new_header = (
                    f"{line.strip()} diff_virial_11 diff_virial_22 "
                    f"diff_virial_33 diff_virial_12 diff_virial_13 diff_virial_23\n"
                )
                lines_buffer.append(new_header)
            else:
                if atom_index < n_atoms:
                    # Append the difference values to the existing line
                    diffs_str = " ".join(str(diffs_lists[d][atom_index]) for d in directions_order)
                    new_line = f"{line.strip()} {diffs_str}\n"
                    lines_buffer.append(new_line)
                    atom_index += 1
                else:
                    print(f"Warning: More lines than data at line {line_num}")

        f_out.writelines(lines_buffer)


def process_files(folder, start, end, delta=0.01):
    """
    Reads original states & deformed states with ± strain
    in 6 directions (xx, yy, zz, xy, xz, yz) . 
    Computes finite difference in potential energies, appends them to a new file.
    """
    # Columns in the original, undeformed file
    original_pe_col, original_vol_col = 8, 9  
    # Columns in the deformed files (xxP, xxN, etc.)
    deformed_pe_col, deformed_vol_col = 1, 2

    directions_map = {
        "xx": "11",
        "yy": "22",
        "zz": "33",
        "xy": "12",
        "xz": "13",
        "yz": "23"
    }
    # Finite difference factor for atomic virial calculations
    half_delta_inv = 1.0 / (2.0 * delta)
    
    for i in range(start, end + 1):
        # 1) Original file name
        ori_filename = os.path.join(folder, f"dump.{i}.init")
        if not os.path.exists(ori_filename):
            print(f"File {ori_filename} does not exist. Skipping.")
            continue

        # 2) Read the original frame
        ori_df = read_lammps_data(ori_filename, original_pe_col, original_vol_col)
        # Finite difference factor for atomic virial stress calculations
        denominator = 1.0 / (2 * ori_df.vol * delta) * evA2Gpa

        # 3) Read strained frames
        diff_virial = {}
        for direction, label in directions_map.items():
            P_file = os.path.join(folder, f"dump.{i}.{direction}P")
            N_file = os.path.join(folder, f"dump.{i}.{direction}N")

            if not os.path.exists(P_file) or not os.path.exists(N_file):
                print(f"Warning: Missing files for direction {direction} at index {i}. Skipping.")
                diff_virial[label] = pd.Series([0.0]*len(ori_df))
                continue

            df_P = read_lammps_data(P_file, deformed_pe_col, deformed_vol_col)
            df_N = read_lammps_data(N_file, deformed_pe_col, deformed_vol_col)

            # Compute finite difference: (df_P - df_N) * (1 / (2 * delta))
            diff_virial[label] = compute_diff_virial(df_P, df_N, half_delta_inv)
            # Replace 'half_delta_inv' with 'denominator' if to compute atomic stress
            # diff_virial[label] = compute_diff_virial(df_P, df_N, denominator)

        # 4) Write out the new file (with appended diff columns)
        output_folder = "dumpVirial_test"
        output_filename = os.path.join(output_folder, f"dump.{i}.cfg")
        write_cfg_with_diff(ori_filename, output_filename, diff_virial)

       #  print(f"Processed and written to {output_filename}")


if __name__ == "__main__":
    # Example usage
    process_files(folder="dumpDiff_test/", start=0, end=50, delta=0.01)
