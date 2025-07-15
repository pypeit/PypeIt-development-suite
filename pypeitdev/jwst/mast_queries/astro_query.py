from astroquery.mast import MastMissionsClass


# Create object and initialize mission to JWST, 
token = "985000d273fd4ed29b24205cea69c7e0"
MastClass = MastMissionsClass(mission='JWST')

# Then login with your token
# go to https://auth.mast.stsci.edu/token to get a new one. 
MastClass.login(token=token)


# Query for datasets
datasets = MastClass.query_criteria(program='9180',
                            opmode='FIXEDSLIT',  # lamp operating mode
                            exp_type='NRS_FIXEDSLIT',  # exposure type
                            productLevel='1b')  # product level (uncalibrated)
 
# Fetch products
products = MastClass.get_unique_product_list(datasets)
 
# Filter for the desired files
# This is how you would limit product files to a certain file type
filtered = MastClass.filter_products(products, file_suffix=['_uncal', '_rate', '_msa'], extension='fits')
#filtered = MastClass.filter_products(products, file_suffix=['_rate', '_msa'], extension='fits')


# Download files
download_dir = '/Users/joe/jwst_redux/Raw/NIRSPEC_FS/9180'
MastClass.download_products(filtered, download_dir=download_dir, verbose=True)