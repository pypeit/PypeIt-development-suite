from astroquery.mast import MastMissionsClass

# Create object and initialize mission to JWST
MastClass = MastMissionsClass(mission='JWST')
 
# Query for datasets
datasets = MastClass.query_criteria(program='3543',
                            opmode='MSASPEC',  # lamp operating mode
                            exp_type='NRS_MSASPEC',  # exposure type
                            productLevel='1b')  # product level (uncalibrated)
 
# Fetch products
products = MastClass.get_unique_product_list(datasets)
 
# Filter for the desired files
# This is how you would limit product files to a certain file type
#filtered = MastClass.filter_products(products, file_suffix=['_uncal', '_rate', '_msa'], extension='fits')
filtered = MastClass.filter_products(products, file_suffix=['_rate', '_msa'], extension='fits')


# Download files
download_dir = '/Users/joe/jwst_redux/Raw/NIRSPEC_MSA/3543'
MastClass.download_products(filtered, download_dir=download_dir, verbose=True)