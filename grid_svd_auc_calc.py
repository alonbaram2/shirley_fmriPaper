import numpy as np
import scipy.io
import os.path


class Svd_AUC_Grid:
    def __init__(self, grid_cells=np.empty(2), within_left_grid=0, between_left_grid=0, permuted_g_auc_dif=np.empty(2), permuted_g_auc_between=np.empty(2),n_permutation=1000, n_grid=41):

        self.grid_cells = grid_cells[:, :n_grid, :]
        self.within_left_grid = within_left_grid
        self.between_left_grid = between_left_grid
        self.permuted_g_auc_dif = permuted_g_auc_dif
        self.permuted_g_auc_between = permuted_g_auc_between
        self.num_permutation = n_permutation

    """
    permuted_g_auc - permute grid cells auc (we permute the cells in the PCs)

    """

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

    def grid_cells_SVDs(self, area_name):

        """

        Parameters:
        area_name (str): whether to calculate the svd analysis over grid/place cells/permuted data
        trials(np.array): if not empty it should be the permuted data

        Returns:
        tuple: A tuple (cum_var_left, cum_var_right), where cum_var_left and cum_var_right are lists of cumulative
               variances for left and right singular vectors for each trial.

        Performs Singular Value Decomposition (SVD) analysis on place grid cell data.
        It then performs the SVD analysis for the first trial, and for
        each subsequent trial, it calculates the cumulative variances using the left and right singular vectors
        obtained from the first trial's SVD. The cumulative variances for each trial are then returned.
        """
        #
        # Note that below I just selected the first 41 cells  to make grid and place cell area under the curve comparable (as this is the total number of grid cells but there are more place cells in the data).
        # This is not ideal, we should do some shuffling instead.
        #

        trials = self.grid_cells

        t_u, t_v = self.svd_analysis(trials[0])  # SVD of the first trial
        print(t_u.shape)
        if area_name == 'permuted':
            t_u = t_u[np.random.permutation(t_u.shape[0])]

        cum_var_left = []
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

    def cal_auc_real(self):
        """
        This function calculate grid and place cells auc over the real data
        """
        between_left_grid, within_left_grid, between_right_grid, within_right_grid = self.grid_cells_SVDs('grid')

        print(
            f"within_left_grid sum: {np.sum(within_left_grid / within_left_grid.shape[0])} and shape: {within_left_grid.shape} \
              between_left_grid sum: {np.sum(between_left_grid / between_left_grid.shape[0])} and shape: {between_left_grid.shape}")
        self.within_left_grid = np.sum(within_left_grid) / within_left_grid.shape[0]
        self.between_left_grid = np.sum(between_left_grid) / within_left_grid.shape[0]
        auc_dif_grid = np.sum(within_left_grid - between_left_grid) / within_left_grid.shape[0]
        print(auc_dif_grid)
        return auc_dif_grid

    def cal_auc_permuted(self):
        """
        This function calculate the auc of the permuted data

        Return: The auc difference
        """
        between_left, within_left, between_right, within_right = self.grid_cells_SVDs(
            'permuted')
        auc_dif = np.sum(within_left - between_left) / within_left.shape[0]
        between_left_auc_dif = self.within_left_grid - np.sum(between_left) / within_left.shape[0]
        return auc_dif, between_left_auc_dif


    def cal_auc_permuted_vec(self):
        """
        This function calculate the dif-auc over the permuted/sampled data
        If name_permutation = 'grid_place' the permutation is over grid and place cells
        if name_permutation = 'place' the place cells are sampled and the difference of place cells
            auc is calculated
        """
        auc_dif_arr = np.array([self.cal_auc_permuted() for x in np.arange(self.num_permutation)])
        print(f"auc_dif_arr shape: {auc_dif_arr.shape}")
        self.permuted_g_auc_dif = auc_dif_arr[:, 0]
        self.permuted_g_auc_between = auc_dif_arr[:, 1]


    def create_permuted_auc_dist(self, name_permutation='between'):
        """
        This function calculates the permuted dif-auc distribution.
        Parameters: name_permutation (str)
            If name_permutation = 'grid_place' the permutation is over grid and place cells
            if name_permutation = 'place' the place cells are sampled and the difference of place cells
                auc is calculated

        Return: the cdf
        """
        if name_permutation == 'between':
            hist, bin_edges = np.histogram(self.permuted_g_auc_between, bins=100)
        elif name_permutation == 'dif':
            hist, bin_edges = np.histogram(self.permuted_g_auc_dif, bins=100)

        x = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        prob = hist / np.sum(hist)
        return x, np.cumsum(prob)



