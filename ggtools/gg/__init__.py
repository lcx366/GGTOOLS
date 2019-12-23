'''
ggtools gg subpackage

This subpackage defines the following functions:

# ====================== fitting function ==================== #

func - Define a linear function f(x) = a0 + a1/T*x to be fitted, where a0 and a1 are parameters for inter

# ====================== estinate lovebums =================== #

lovebums - Estimate Load Love Numbers(LLNs)

# ===================== plotting functions =================== #

plot_at_northpole - Plot grid data at the northpole.  

# ===================== mascon functions ===================== #

mascon_download - Download GRACE mascon data from https://neptune.gsfc.nasa.gov/uploads/grace/mascons_2.4/

read_mascon - Read the GRACE mascon files.

# ========== Compare DDK filter and Gaussian filter ========== #

ddk_gausian - Given a specific type of DDK filter and the maximum SHC degree number, evaluate the 'equivalent' Gaussian filter radius. 

# ================= Tikhonov Regularization ================== #

solve_lambda_x - Minimize the Lagrangian L(x,lambda) = || A*x - y||**2 + lambda*||x||**2 based on the idea of Tikhonov Regularization. Estimate the parameters x and Lagrange multiplier roughly.

L_curve - Minimize the Lagrangian L(x,lambda) = || A*x - y||**2 + lambda*||x||**2 based on the idea of Tikhonov Regularization. Estimate the parameters x and Lagrange multiplier accurately.

# ===================== Filter functions ===================== #

filter_ddk - DDK filter used to attenuate noise described as striping patterns in GRACE GSM data. 

filter_gaussian - Gaussian filter used to attenuate noise described as striping patterns in GRACE GSM data.

filter_gaussian_inverse - Inversion of Gaussian filter. It is used to recover the signal in the process of leakage correction in GRACE GSM data.

# ================= Signal leakage correction ================ #

scale_factor - Estimate the scale factor(gain factor) used in signal leakage correction. 

forward_model_initial - Expand the grid data to spherical harmonic coefficients and perform Gaussian filtering, then transfer it back to the grid data. 

forward_model - Iterative Forward Modeling used to perform signal leakage correction in GRACE data processing.

space_domain - Space domain method used to perform signal leakage correction in GRACE data processing.

spectral_domain - Spectrum domain method used to perform signal leakage correction in GRACE data processing.

# ==================== Least square method =================== #

lsqm - Linearly fit f(x) = a0 + a1/T*x using the Least Square algorithm.

wlsqm - Linearly fit f(x) = a0 + a1/T*x using the Weighted Least Square algorithm.

ilsqm - Linearly fit f(x) = a0 + a1/T*x using the Iterative Least Square method. The 3-sigma rule is used to eliminate outliers.

iwlsqm - Linearly fit f(x) = a0 + a1/T*x using the Iterative Weighted Least Square algorithm. The 3-sigma rule is used to eliminate outliers.

# ==================== GRACE data utilitues ================== #

parse_gsm_filename - Parse GRACE GSM filenames.

print_gsm_date_coverage - Print the date coverage for the GRACE GSM data from isdcftp.gfz-potsdam.de

gsm_download - Download GRACE GSM data from isdcftp.gfz-potsdam.de

parse_gsm_file - Parse the GRACE GSM Level-2 file.

read_gsm - Read the GRACE GSM Level-2 files.

gsm_average - Combine the (deaveraged) GSM solution from multiple institutions into an integrated one.

# ==================== GLDAS data utilitues ================== #

gldas_download - Download the GLDAS grid data over a period defined by the start date and end date and its documentation from urs.earthdata.nasa.gov

read_gldas - Read the GLDAS files into a GLDAS class instance. 

regular_gldas - Normalize the GLDAS grid data to meet the requirements of spherical harmonic expansion with pyshtools based on the sampling theorem of Driscoll and Healy (1994).

lsm - Calculate some land surface quantities from GLDAS grid data, such as Terrestrial Water Storage Changes(TWSC).

landmask - Establish a land window function based on the global terrain data ETOPO5.

# =================== SLR C20 data utilitues ================= #

slr_c20_download - Download SLR C20 data from isdcftp.gfz-potsdam.de

read_slr_c20 - Read the SLR C20 file.

# ========================= Utilitues ======================== #

print_error_grace - If source, D, and RL do not meet the input requirements, an error message will be printed.

print_error_gldas - If source and res do not meet the input requirements, an error message will be printed.

month2int - Given a list of month list, translate it to an array of month sequence.

med(x) - Calculate the middle divisor.

solid_angle_ratio - Calculate the ratio of a solid angle of an ellipse on a sphere to 4pi.

crop_region - Crop the global grid data to a region of interested.

yx2latlon - Transform index to latitudes and longitudes.

latlon2yx - Transform latitudes and longitudes to index.

generate_mascons - Create mascons within a study area using a set of points or polygons.
'''  
from .gsm_utils import print_gsm_date_coverage,gsm_download,read_gsm,gsm_average
from .static_models import static_download
from .slr_utils import slr_c20_download,read_slr_c20
from .gldas_utils import gldas_download,read_gldas
from .plot import plot_at_northpole
from .ddk_gaussian import ddk_gaussian
from .leakage import scale_factor
from .lcurve import L_curve
from .utils import generate_nodes,solid_angle_ratio
from .mascon import mascon_download,read_mascon
from .ddk_gaussian import ddk_gaussian
from .landmask import landmask