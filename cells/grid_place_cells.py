import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
from scipy import interpolate
from scipy.ndimage import gaussian_filter

DATA_PATH = '/Users/veronikasamborska/Desktop/grid_and_place_cells/'

def interpolate_and_filter_cells(trials):
    """
    Interpolates missing values and applies a Gaussian filter to a list of trials. 
    Note that there are five trials, 1st and 5th are trials when animals are in the same arena, 2nd, 3rd and 4th are from a different arena.
    

    Parameters:
    trials (array-like): A list of 2D numpy arrays where each array represents a trial. 
        The values in the array represent some measurement (e.g., activity of a neuron) at each point.
        Missing values are represented as NaNs.

    Returns:
    numpy.ndarray: A list of flattened, filtered trials in which the original NaN values have been interpolated 
        and a Gaussian filter has been applied. Each row in the returned array represents a flattened, 
        filtered trial.

    For each trial, the function operates as follows:
    1. It identifies valid (i.e., non-NaN) values and their coordinates.
    2. It uses these values and coordinates to create an interpolator function using scipy's LinearNDInterpolator. 
        This function can estimate the missing (i.e., NaN) values based on the valid ones.
    3. It applies the interpolator to the original image, filling in NaN values. The result is a 2D array of the 
        same shape as the input with no NaN values.
    4. It applies a Gaussian filter to the interpolated image to smooth it.
    5. It flattens the 2D filtered image array to a 1D array.
    6. The final flattened, filtered image array is then appended to the list of processed trials.

    The function processes each trial in the same way and returns the list of processed trials as a 2D numpy array.
    
    """
    interpolated_trials = []
    for trial in trials:
        image = trial
        valid_mask = ~np.isnan(image)
        coords = np.array(np.nonzero(valid_mask)).T
        values = image[valid_mask]
        it = interpolate.LinearNDInterpolator(coords, values, fill_value=0)
        filled = it(list(np.ndindex(image.shape))).reshape(image.shape)
        filled = gaussian_filter(filled, 2)
        interpolated_trials.append(filled.flatten())
    return np.asarray(interpolated_trials)


def svd_and_manipulation(data, mean):
    """
    Transposes and centers a given data set around a given mean. 
    This is done to account for different baseline firing rates in EC and HPC.  

    Parameters:
    data (numpy.ndarray): The input data to be transposed and centered. This should be a 2D numpy array where 
                          each row represents a trial and each column represents a feature.
    mean (float or numpy.ndarray): The mean value(s) to be subtracted from the data. If mean is a single value, 
                                   it's subtracted from all elements in the data. If mean is an array, it should 
                                   be the same shape as data after transposition, so that element-wise subtraction can occur.

    Returns:
    numpy.ndarray: The manipulated data that has been transposed and mean-centered. It maintains the same dimensions 
                   as the original data.

    The function operates by first transposing the input data, so that each column represents a trial and each row represents 
    a feature. It then subtracts the given mean from this transposed data, centering the data. The function finally returns 
    the transposed, centered data.
    """
    data = np.transpose(data) - mean
    return np.transpose(data)


def svd_analysis(data):
    """
    Performs Singular Value Decomposition (SVD) analysis on the given data.

    Parameters:
    data (numpy.ndarray): The input data to perform SVD on. This should be a 2D numpy array where 
                          each row represents a trial and each column represents a feature.

    Returns:
    tuple: A tuple (t_u, t_v) where t_u and t_v are the transposed left singular vectors and right singular vectors, respectively.

    This function performs SVD on the input data using the numpy's linalg.svd method with full_matrices set to True, 
    which means it produces u and vh matrices with shapes that are compatible with original data shape. 
    The left singular vectors (u) and the right singular vectors (vh) are then transposed and returned.
    """
    u, s, vh = np.linalg.svd(data, full_matrices = True)
    t_u = np.transpose(u)
    t_v = np.transpose(vh)
    return t_u, t_v


def calculate_cumulative_variances(data, t_u, t_v):
    """
    Calculates the cumulative variances for the transformed input data.

    Parameters:
    data (numpy.ndarray): The input data for which the cumulative variances are to be calculated.
    t_u (numpy.ndarray): The transposed left singular vectors obtained from SVD.
    t_v (numpy.ndarray): The transposed right singular vectors obtained from SVD.

    Returns:
    tuple: A tuple (cum_var_x, cum_var_y) where cum_var_x and cum_var_y are the cumulative variances for 
           the left and right singular vectors transformations, respectively.

    The function first performs a dot product between the left and right singular vectors and the input data, 
    calculating the sum of squares along the axis 1 (rows). It then calculates the cumulative sum of the variance, 
    normalized by the square root of the data's shape[0] (number of rows), and finally normalizes the cumulative 
    variance by its last value, to provide a percentage-like format.
    """
    x = np.linalg.multi_dot([t_u, data])
    var_x = np.sum(x**2, axis = 1)
    cum_var_x = np.cumsum(var_x)/np.sqrt(data.shape[0])
    cum_var_x = cum_var_x/cum_var_x[-1]

    y = np.linalg.multi_dot([data, t_v])
    var_y = np.sum(y**2, axis = 1)
    cum_var_y = np.cumsum(var_y)/np.sqrt(data.shape[0])
    cum_var_y = cum_var_y/cum_var_y[-1]

    return cum_var_x, cum_var_y


def place_grid_cells_SVDs(file, dictionary, color='green'):
    
    """
    Performs Singular Value Decomposition (SVD) analysis on place grid cell data.

    Parameters:
    file (str): The filename containing the place grid cell data.
    dictionary (str): The dictionary key to retrieve the place grid cell data from the .mat file.
    color (str, optional): The color to use for plotting. Default is 'green'.

    Returns:
    tuple: A tuple (cum_var_left, cum_var_right), where cum_var_left and cum_var_right are lists of cumulative 
           variances for left and right singular vectors for each trial.

    The function loads the .mat file, and for each trial, it first interpolates and filters the place cells, 
    then applies SVD and manipulates the data. It then performs the SVD analysis for the first trial, and for 
    each subsequent trial, it calculates the cumulative variances using the left and right singular vectors 
    obtained from the first trial's SVD. The cumulative variances for each trial are then returned.
    """
    load_cells = scipy.io.loadmat(DATA_PATH + file)
    
    #
    # Note that below I just selected the first 41 cells  to make grid and place cell area under the curve comparable (as this is the total number of grid cells but there are more place cells in the data). 
    # This is not ideal, we should do some shuffling instead. 
    #
    
    extract_cells = load_cells[dictionary][:41]  

    trials = [interpolate_and_filter_cells(extract_cells[:, i]) for i in range(5)]

    mean = np.mean(np.concatenate(trials, axis=1), axis=1)
    trials = [svd_and_manipulation(t, mean) for t in trials]

    t_u, t_v = svd_analysis(trials[0]) # SVD of the first trial 
    
    cum_var_left = [];  cum_var_right = []
    for i, t in enumerate(trials[1:]):  # SVD trials 2,3,4,5 
        cum_var_x, cum_var_y = calculate_cumulative_variances(t, t_u, t_v)
        cum_var_left.append(cum_var_x)
        cum_var_right.append(cum_var_y)
    
    between_left = np.mean(cum_var_left[:3],0) # 2,3 and 4 is a different arena hence this is between
    within_left = cum_var_left[-1] # Last trial is in the same arean as the first hence here this is within
    
    between_right = np.mean(cum_var_right[:3],0) # 2,3 and 4 is a different arena hence this is between
    within_right = cum_var_right[-1] # Last trial is in the same arean as the first hence here this is within

    return between_left, within_left, between_right, within_right
  
    

def plot_cumulative_variances(cum_var_within, cum_var_between, sub = 1, color='green', area = 'Grid cells', title = 'Cells'):
    """Plot cumulative variances."""

    plt.subplot(1,2, sub)
    plt.plot(cum_var_within, color=color, alpha=0.7, label = "Within arenas " + area)
    plt.plot(cum_var_between, color=color, alpha=0.7, linestyle = '--', label = "Between arenas -" + area)
    plt.legend()

    plt.xlabel('# of Singular Vectors')
    plt.ylabel('% Variance Explained')
    sns.despine()
    plt.title(title)
 


def run():
    """This function just runs everything and plots within and between variance explained"""
    
    data_files = [
        {'file': 'smthRm_place_all_animals.mat', 'dictionary': 'smthRm_place_all_animals', 'color': 'green'},
        {'file': 'smthRm_grid_all_animals.mat', 'dictionary': 'smthRm_grid_all_animals', 'color': 'black'}
    ]
    
    between_left_PC, within_left_PC, between_right_PC, within_right_PC = place_grid_cells_SVDs(**data_files[0])
    between_left_GR, within_left_GR, between_right_GR, within_right_GR = place_grid_cells_SVDs(**data_files[1])
    
    plt.figure(figsize=(12, 6))

    plot_cumulative_variances(within_left_GR, between_left_GR, sub = 1, color='black', area  = "grid cells", title = "Cells")
    plot_cumulative_variances(within_left_PC, between_left_PC, sub = 1, color='green', area = "place cells",title = "Cells")

        
    plot_cumulative_variances(within_right_GR, between_right_GR, sub = 2, color='black', area  = "grid cells", title = "Location bins")
    plot_cumulative_variances(within_right_PC, between_right_PC, sub = 2, color='green', area  = "place cells", title = "Location bins")

        
if __name__ == '__main__':
    run()