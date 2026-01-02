import cdsapi
import xarray as xr
import pandas as pd
import os

# --- CONFIGURATION TEMPLATE ---
# If you didn't set up the .cdsapirc file, fill these in:
CDS_URL = 'https://cds.climate.copernicus.eu/api'
CDS_KEY = '4c15ef70-ac4a-472b-9313-22fa91387ebb'  # e.g., '12345:88888888-4444-....'
# ------------------------------

def check_file_validity(filepath):
    """
    Checks if the file is a valid NetCDF by trying to read the header.
    Often, failed downloads save as text files containing error HTML.
    """
    try:
        # 1. Check size (arbitrary small threshold like 1KB)
        size_bytes = os.path.getsize(filepath)
        if size_bytes < 1000:
            with open(filepath, 'r') as f:
                content = f.read()
            print(f"\n[!] File is too small ({size_bytes} bytes). Content preview:\n{content[:200]}")
            return False
            
        # 2. Try opening with netCDF4 library directly (lighter than xarray)
        import netCDF4
        with netCDF4.Dataset(filepath) as nc:
            pass # Just opening it is enough to test validity
        return True
        
    except Exception as e:
        print(f"\n[!] File check failed: {e}")
        return False

def download_with_retry(client, dataset, request, target_file, max_retries=3):
    attempt = 1
    while attempt <= max_retries:
        try:
            print(f"Attempt {attempt}/{max_retries}...")
            client.retrieve(dataset, request, target_file)
            
            # Immediate validity check
            if check_file_validity(target_file):
                print("Download successful and verified.")
                return True
            else:
                print("Download 'finished' but file is invalid (likely server error). Retrying...")
                
        except Exception as e:
            print(f"Network/API Error on attempt {attempt}: {e}")
            
        # Exponential backoff
        wait_time = 15 * attempt
        print(f"Waiting {wait_time}s...")
        time.sleep(wait_time)
        attempt += 1
    
    raise RuntimeError("Failed to download valid data after multiple attempts.")

def fetch_iowa_core_rain(year, month, output_file='iowa_core_precip.nc'):
    
    if 'YOUR_UID' in CDS_KEY:
         c = cdsapi.Client() # Relies on .cdsapirc
    else:
        c = cdsapi.Client(url=CDS_URL, key=CDS_KEY)

    # --- OPTIMIZED BOUNDING BOX: CENTRAL IOWA ---
    # ~600 grid points vs ~13,000 previously
    # Covers: Des Moines, Ames, Fort Dodge
    iowa_core_bbox = [43.5, -95.0, 41.5, -92.0] # N, W, S, E

    if isinstance(month, str): month = [month]

    request_params = {
        'variable': 'total_precipitation',
        'year': year,
        'month': month,
        'day': [str(d).zfill(2) for d in range(1, 32)],
        'time': [f"{h:02d}:00" for h in range(24)],
        'area': iowa_core_bbox,
        'format': 'netcdf',
        # --- KEY CHANGE HERE ---
        # Resample to 0.25 degrees (approx 27km)
        # This reduces data volume by ~84%
        'grid': [0.25, 0.25], 
        # -----------------------
    }

    print(f"Requesting data for {year}-{month} (Central Iowa Patch)...")
    download_with_retry(c, 'reanalysis-era5-land', request_params, output_file)

    print("Opening dataset...")
    # Explicitly using netcdf4 engine to avoid ambiguity
    ds = xr.open_dataset(output_file, engine='netcdf4')
    
    # Convert to DataFrame
    df = ds.to_dataframe().reset_index()
    
    # Clean up: Convert meters to mm for easier reading
    df['tp_mm'] = df['tp'] * 1000
    
    ds.close()
    print(f"Success! DataFrame shape: {df.shape}")
    return df

if __name__ == "__main__":
    # Test Run
    try:
        df = fetch_iowa_core_rain('2012', '07')
        print(df.head())
    except Exception as e:
        print(f"Fatal Error: {e}")