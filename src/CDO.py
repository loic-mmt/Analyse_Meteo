from cdo import Cdo
import os

def compute_general_climatology(
    input_file, 
    weights_file, 
    output_path, 
    year_range=None, 
    mode="yearly", 
    selected_months=None, 
    selected_days=None, 
    selected_hours=None,
    variable_name="t2m"
):
    """
    Args:
        input_file (str): Path to the global .nc file.
        weights_file (str): Path to the weights/mask .nc file.
        output_path (str): Where to save the final result.
        year_range (range/list): e.g., range(1950, 2026).
        mode (str): 'yearly', 'monthly', 'daily', or 'total'.
        selected_months (list): e.g., [6, 7, 8].
        selected_days (list): e.g., [1, 15] or None.
        selected_hours (list): e.g., [18] or range(0, 24).
        variable_name (str): Variable to process (default 't2m').
    """
    cdo = Cdo()
    
    
    # Start with the input file
    ops = f"{input_file}"
    
    # A. Select Variable
    ops = f"-selname,{variable_name} {ops}"
    
    # B. Select Years
    if year_range:
        y_min, y_max = min(year_range), max(year_range)
        ops = f"-selyear,{y_min}/{y_max} {ops}"
        
    # C. Select Months
    if selected_months:
        # Convert list [6, 7] to string "6,7"
        m_str = ",".join(map(str, selected_months))
        ops = f"-selmon,{m_str} {ops}"
        
    # D. Select Days (Optional)
    if selected_days:
        d_str = ",".join(map(str, selected_days))
        ops = f"-selday,{d_str} {ops}"

    # E. Select Hours (Optional)
    if selected_hours is not None:
        # Handle range objects or lists
        h_str = ",".join(map(str, list(selected_hours)))
        ops = f"-selhour,{h_str} {ops}"

    # Convert Kelvin to Celsius (Data - 273.15) BEFORE averaging
    ops = f"-subc,273.15 {ops}"

    # Map the 'mode' symbol to the correct CDO operator
    if mode == "yearly":
        ops = f"-yearmean {ops}"
    elif mode == "monthly":
        ops = f"-monmean {ops}"
    elif mode == "daily":
        ops = f"-daymean {ops}"
    elif mode == "total":
        ops = f"-timmean {ops}"
    else:
        print(f"Unknown mode '{mode}'.")

    # We prepare the weight file: treating 0.0 as Missing (NaN)
    weight_ops = f"-setctomiss,0 -gtc,0 -selname,final_weights {weights_file}"
    
    # Final Command: Multiply (Data Chain) * (Weight Chain)
    try:
        cdo.mul(input=f"{ops} {weight_ops}", output=output_path)
        print(f"Saved to {output_path}")
        
    except Exception as e:
        print(e)

if __name__ == "__main__":
    
    # Define your paths
    INPUT_FILE = "src/era5_global_1950_2025_basic.nc"
    WEIGHTS_FILE = "data/masks/weights_prop_basic.nc"
    
    # 1. Standard Yearly Climatology (All year)
    compute_general_climatology(
        INPUT_FILE, WEIGHTS_FILE, 
        output_path="output/climatology_yearly.nc",
        year_range=range(1950, 2026),
        mode="yearly"
    )