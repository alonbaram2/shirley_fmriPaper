import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os.path
from scipy import interpolate
from scipy.ndimage import gaussian_filter
import os

from grid_svd_auc_calc import Svd_AUC_Grid
DATA_PATH = 'cells/grid_and_place_cells/grid_and_place_cells'


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


def d_mean_manipulation(data, mean):
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


def data_preprocessing(file, dictionary, color='green'):
    """
    Performs Singular Value Decomposition (SVD) analysis on place grid cell data.

    Parameters:
    file (str): The filename containing the place grid cell data.
    dictionary (str): The dictionary key to retrieve the place grid cell data from the .mat file.
    color (str, optional): The color to use for plotting. Default is 'green'.

    Returns:
    numpy.ndarray: The manipulated data that has been transposed and mean-centered. It maintains the same dimensions
                   as the original data.

    The function loads the .mat file, and for each trial, it interpolates and filters the place cells.
    """
    print(os.getcwd())
    print(os.path.join(DATA_PATH, file))
    load_cells = scipy.io.loadmat(os.path.join(os.getcwd(), DATA_PATH, file))

    #
    # Note that below I just selected the first 41 cells  to make grid and place cell area under the curve comparable (as this is the total number of grid cells but there are more place cells in the data).
    # This is not ideal, we should do some shuffling instead.
    #

    extract_cells = load_cells[dictionary][:41]
    print(extract_cells.shape)
    trials = [interpolate_and_filter_cells(extract_cells[:, i]) for i in range(5)]

    mean = np.mean(np.concatenate(trials, axis=1), axis=1)
    trials = [d_mean_manipulation(t, mean) for t in trials]

    return np.array(trials)


def run():
    """This function runs the preprcessing"""

    data_files = [
        {'file': 'smthRm_place_all_animals.mat', 'dictionary': 'smthRm_place_all_animals', 'color': 'green'},
        {'file': 'smthRm_grid_all_animals.mat', 'dictionary': 'smthRm_grid_all_animals', 'color': 'black'}
    ]

    place_cells = data_preprocessing(**data_files[0])
    grid_cells = data_preprocessing(**data_files[1])

    return place_cells, grid_cells

def my_func(a):
    print(a.shape)

if __name__ == '__main__':
    place_cells, grid_cells = run()
    print(f"grid cells array shape {grid_cells.shape}")
    num_permutation = 1000
    permuted_auc_dif_file_name = 'cells/grid_permute_auc_dif' + str(num_permutation)
    permuted_auc_bet_file_name = 'cells/grid_permute_auc_between' + str(num_permutation)
    if os.path.exists(permuted_auc_dif_file_name + '.npy'):
        print('auc permutation file exists')
        permuted_g_auc_dif = np.load(permuted_auc_dif_file_name + '.npy')
        permuted_g_auc_between = np.load(permuted_auc_bet_file_name + '.npy')
        Svd_AUC_obj = Svd_AUC_Grid(grid_cells=grid_cells, permuted_g_auc_dif=permuted_g_auc_dif, permuted_g_auc_between=permuted_g_auc_between)
    else:
        Svd_AUC_obj = Svd_AUC_Grid(grid_cells=grid_cells, n_permutation=num_permutation)
        Svd_AUC_obj.cal_auc_permuted_vec()
        np.save(permuted_auc_dif_file_name, Svd_AUC_obj.permuted_g_auc_dif)
        np.save(permuted_auc_bet_file_name, Svd_AUC_obj.permuted_g_auc_between)
    Svd_AUC_obj.cal_auc_real()
    print(Svd_AUC_obj.within_left_grid, Svd_AUC_obj.between_left_grid)
    print(Svd_AUC_obj.permuted_g_auc_dif.shape)
    x_dif, cdf_dif = Svd_AUC_obj.create_permuted_auc_dist(name_permutation='dif', n_bin=250)
    x_between, cdf_between = Svd_AUC_obj.create_permuted_auc_dist(name_permutation='between', n_bin=250)
    x_between_dif, cdf_between_dif = Svd_AUC_obj.create_permuted_auc_dist(name_permutation='between_dif', n_bin=200)
    grid_auc_effect = Svd_AUC_obj.within_left_grid - Svd_AUC_obj.between_left_grid
    n_auc_effect = np.sum(x_between_dif < grid_auc_effect)
    p_val_auc_grid_dif = cdf_between_dif[n_auc_effect]
    print(f"p_values: {p_val_auc_grid_dif}")
    plt.figure(1)
    plt.plot(x_dif, cdf_dif)
    plt.figure(2)
    plt.plot(x_between, cdf_between)
    plt.figure(3)
    plt.plot(x_between_dif, cdf_between_dif)
    plt.figure(4)
    plt.plot(Svd_AUC_obj.within_left_grid_plot, 'k', label='same arena')
    plt.plot(Svd_AUC_obj.between_left_grid_plot, 'k-.', label='different arena')
    plt.plot(np.arange(0, 1, 1/41), 'k:', label='random')
    plt.legend()
    plt.show()




