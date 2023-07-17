import numpy as np
import scipy.io
import os.path
from scipy import interpolate
from scipy.ndimage import gaussian_filter

class Svd_AUC:
    def __init__(self, grid_cells=np.empty(2), place_cells=np.empty(2), num_permutations=10):
        self.grid_cells = grid_cells[:, :40, :]
        self.place_cells = place_cells
        self.num_per = num_permutations
        self.permuted_data = np.empty((self.num_per,) + self.grid_cells.shape)
    def svd_analysis(self, data):
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
        u, s, vh = np.linalg.svd(data, full_matrices=True)
        t_u = np.transpose(u)
        t_v = np.transpose(vh)
        return t_u, t_v

    def calculate_cumulative_variances(self, data, t_u, t_v):
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
        var_x = np.sum(x ** 2, axis=1)
        cum_var_x = np.cumsum(var_x) / np.sqrt(data.shape[0])
        cum_var_x = cum_var_x / cum_var_x[-1]

        y = np.linalg.multi_dot([data, t_v])
        var_y = np.sum(y ** 2, axis=1)
        cum_var_y = np.cumsum(var_y) / np.sqrt(data.shape[0])
        cum_var_y = cum_var_y / cum_var_y[-1]

        return cum_var_x, cum_var_y

    def place_grid_cells_SVDs(self, area_name):

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
        print(os.getcwd())

        #
        # Note that below I just selected the first 41 cells  to make grid and place cell area under the curve comparable (as this is the total number of grid cells but there are more place cells in the data).
        # This is not ideal, we should do some shuffling instead.
        #
        if area_name=='grid':
            trials = self.grid_cells
        if area_name=='place':
            trials = self.place_cells

        t_u, t_v = self.svd_analysis(trials[0])  # SVD of the first trial

        cum_var_left = [];
        cum_var_right = []
        for i, t in enumerate(trials[1:]):  # SVD trials 2,3,4,5
            cum_var_x, cum_var_y = self.calculate_cumulative_variances(t, t_u, t_v)
            cum_var_left.append(cum_var_x)
            cum_var_right.append(cum_var_y)

        between_left = np.mean(cum_var_left[:3], 0)  # 2,3 and 4 is a different arena hence this is between
        within_left = cum_var_left[-1]  # Last trial is in the same arean as the first hence here this is within

        between_right = np.mean(cum_var_right[:3], 0)  # 2,3 and 4 is a different arena hence this is between
        within_right = cum_var_right[-1]  # Last trial is in the same arean as the first hence here this is within

        return between_left, within_left, between_right, within_right


    def cal_auc(self):
        between_left_grid, within_left_grid, between_right_grid, within_right_grid = self.place_grid_cells_SVDs('grid')
        between_left_place, within_left_place, between_right_grid, within_right_place= self.place_grid_cells_SVDs('place')
        auc_dif_grid = np.sum(within_left_grid - between_left_grid)
        auc_dif_place = np.sum(within_left_place - between_left_place)
        print(auc_dif_grid, auc_dif_place)

    def permute_data(self, num_samples=10):
        array1 = self.place_cells
        array2 = self.grid_cells
        num_dims1 = array1.shape[1]
        num_dims2 = array2.shape[1]
        print(num_dims1, num_dims2)
        half_dims = min(num_dims1, num_dims2) // 2

        sampled_arr = []
        for _ in range(num_samples):
            # Sample half of the dimensions from each array
            sampled_indices1 = np.random.choice(num_dims1, size=half_dims, replace=False)
            sampled_indices2 = np.random.choice(num_dims2, size=half_dims, replace=False)

            # Select the sampled dimensions from each array
            sampled_array1 = array1[:, sampled_indices1, :]
            sampled_array2 = array2[:, sampled_indices2, :]

            # Concatenate the sampled dimensions from both arrays
            sampled_arr.append(np.concatenate((sampled_array1, sampled_array2), axis=1))
            self.permuted_data = np.array(sampled_arr)


