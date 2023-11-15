import numpy as np
import scipy.io
import os.path


class Svd_AUC:
    def __init__(self, grid_cells=np.empty(2), place_cells=np.empty(2), permuted_gp_auc=np.empty(2), sampled_p_data=np.empty(2), sampled_p_auc=np.empty(2), n_grid=41):
        self.num_real_place_cells = n_grid
        self.grid_cells = grid_cells[:, :n_grid, :]
        self.place_cells = place_cells
        self.permuted_gp_data = np.empty(self.grid_cells.shape)
        self.sampled_p_data = sampled_p_data
        self.permuted_gp_auc = permuted_gp_auc
        self.sampled_p_auc = sampled_p_auc
    """
    permuted_gp_data - permute grid and place cells
    permuted_p_data - sample place cells
    
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

    def place_grid_cells_SVDs(self, area_name, trials=np.empty(2)):

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
        if area_name == 'grid':
            trials = self.grid_cells
        if area_name == 'place':
            trials = self.place_cells[:self.num_real_place_cells]

        t_u, t_v = self.svd_analysis(trials[0])  # SVD of the first trial

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
        between_left_grid, within_left_grid, between_right_grid, within_right_grid = self.place_grid_cells_SVDs('grid')
        between_left_place, within_left_place, between_right_grid, within_right_place = self.place_grid_cells_SVDs('place')
        print(f"within_left_grid sum: {np.sum(within_left_grid / within_left_grid.shape[0])} and shape: {within_left_grid.shape}")
        auc_dif_grid = np.sum(within_left_grid - between_left_grid) / within_left_grid.shape[0]
        auc_dif_place = np.sum(within_left_place - between_left_place) / within_left_place.shape[0]
        print(auc_dif_grid, auc_dif_place)
        return auc_dif_grid, auc_dif_place

    def cal_auc_permuted(self, permuted_data):
        """
        This function calculate the auc of the permuted data

        Parameters:
            permuted_data (np.array): the permuted data to use for the dif-auc calculation

        Return: The auc difference
        """
        between_left, within_left, between_right, within_right = self.place_grid_cells_SVDs(
            'permuted', permuted_data)
        auc_dif = np.sum(within_left - between_left) / within_left.shape[0]
        return auc_dif

    def cal_auc_permuted_vec(self, name_permutation='grid_place'):
        """
        This function calculate the dif-auc over the permuted/sampled data
        If name_permutation = 'grid_place' the permutation is over grid and place cells
        if name_permutation = 'place' the place cells are sampled and the difference of place cells
            auc is calculated
        """
        if name_permutation == 'grid_place':
            auc_dif_list = [self.cal_auc_permuted(np.squeeze(self.permuted_gp_data[x, :, :, :])) for x in np.arange(self.permuted_gp_data.shape[0])]
            self.permuted_gp_auc = np.array(auc_dif_list)
        if name_permutation == 'place':
            auc_dif_list = [self.cal_auc_permuted(np.squeeze(self.sampled_p_data[x, :, :, :])) for x in np.arange(self.sampled_p_data.shape[0])]
            self.sampled_p_auc = np.array(auc_dif_list)

    def create_permuted_auc_dist(self, name_permutation='grid_place'):
        """
        This function calculates the permuted dif-auc distribution.
        Parameters: name_permutation (str)
            If name_permutation = 'grid_place' the permutation is over grid and place cells
            if name_permutation = 'place' the place cells are sampled and the difference of place cells
                auc is calculated

        Return: the cdf
        """
        if name_permutation == 'grid_place':
            hist, bin_edges = np.histogram(self.permuted_gp_auc, bins=100)
        if name_permutation == 'place':
            hist, bin_edges = np.histogram(self.sampled_p_auc, bins=100)
        x = 0.5*(bin_edges[1:] + bin_edges[:-1])
        prob = hist / np.sum(hist)
        return x, np.cumsum(prob)

    def permute_grid_place_data(self, num_samples=10):

        """
        Creates permuted grid and place cells array
        Parameters: num_samples(int)
        """
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
            self.permuted_gp_data = np.array(sampled_arr)

    def sample_place_cells(self, n_cells=41, num_samples=10):
        """
        This function sampled cells from the place cells data.

        Parameter:
        n_cells: the number of cells to sample each time (int)
        num_samples: the number of samples to create (int)
        """
        arr = self.place_cells
        result = []
        for _ in range(num_samples):
            samples = np.random.choice(arr.shape[1], size=n_cells, replace=False)
            result.append(arr[:, samples, :])
        self.sampled_p_data = np.array(result)
        print(f"sampled place cells data shape: {self.sampled_p_data.shape}")


