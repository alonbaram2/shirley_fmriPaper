import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os.path
from scipy import stats

DATA_PATH = 'voxels_data'
class Svd_AUC_voxels:
    def __init__(self, voxels_data=np.empty((28, 100, 4, 40))):
        """
        voxels_data: array of blob voxels with dimensions [num_subjects  x num_voxels x runs x num_piles]
        """
        self.voxels_data = voxels_data
        self.voxels_loo_average = np.random.rand(28, 100, 4, 40)
        self.subject_svd_curve = np.empty((28, 4, 2, 100))  # num_subject x map's PC x map variance explained only for hex
        self.subject_auc = np.empty((28, 4, 2)) #num_subject x map's PC x map variance explained only for hex

    def cal_loo_average(self):
        """
        Leave the current run out and calculate average betas over all other three runs
        """
        runs_vector = np.arange(4)
        for r in np.arange(4):
            mean_indexes = ~np.isin(runs_vector, np.array([r]))
            print(f"current indexes to average over: {mean_indexes}")
            self.voxels_loo_average[:, :, r, :] = np.mean(self.voxels_data[:, :, mean_indexes, :], axis=2)

    def single_subject_map_eigenvector(self, subject_number=0, sub_voxels_data=np.empty((100, 4, 40)), sub_voxels_loo=np.empty((4, 100, 40))):

        runs_vector = np.arange(4)
        for r in np.arange(4): #go over the runs
            mean_indexes = np.isin(runs_vector, np.array([r]))
            for m_u in np.arange(4):# run over the maps
                data_mat = np.squeeze(sub_voxels_loo[:, mean_indexes, m_u*10:(m_u+1)*10])
                print(f"data_mat shape is: {data_mat.shape}")
                t_u, t_v = self.svd_analysis(data_mat)
                if r == 0 and m_u == 0 and subject_number == 0:
                    print(f"u shape is: {t_u.shape}")
                for m_var in np.arange(2):# only explain hex not cluster
                    cum_var_x, cum_var_y = self.calculate_cumulative_variances(np.squeeze(sub_voxels_data[:, mean_indexes, m_var*10:(m_var+1)*10]), t_u, t_v)
                    print(f"cum_var_x shape: {cum_var_x.shape}")
                    print(f"subject_number, m_u, m_var {subject_number, m_u, m_var}")
                    self.subject_svd_curve[subject_number, m_u, m_var, :] = cum_var_x
                    self.subject_auc[subject_number, m_u, m_var] = np.sum(cum_var_x) / cum_var_x.shape[0]
                if r == 0 and m_u == 0 and subject_number == 0:
                    print(f"u shape is: {t_u.shape}")


    def run_subjects_auc_calculation(self):
        for sub in np.arange(28):
            self.single_subject_map_eigenvector(sub, sub_voxels_data=np.squeeze(self.voxels_data[sub, :, :, :]),
                                                sub_voxels_loo=np.squeeze(self.voxels_loo_average[sub, :, :, :]))
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

    def plot_auc_curve(self):
        averaged_projection = np.squeeze(np.mean(self.subject_svd_curve, axis=0))
        within_map_projection = 0.5 * (averaged_projection[0, 0, :] + averaged_projection[1, 1, :])
        grid_projection = 0.5 * (averaged_projection[0, 1, :] + averaged_projection[1, 0, :])
        same_structure = 0.5*(within_map_projection + grid_projection)
        diff_structure_projection = 0.25 * (averaged_projection[3, 1, :] + averaged_projection[3, 0, :] + \
                                            averaged_projection[2, 1, :] + averaged_projection[2, 0, :])
        plt.plot(within_map_projection)
        plt.plot(grid_projection, ':')
        plt.plot(diff_structure_projection, 'r')
        plt.show()

    def cal_main_effect(self):
        same_map = 0.5*(self.subject_auc[:, 0, 0] + self.subject_auc[:, 1, 1])
        same_structure_dif_map = 0.5*(self.subject_auc[:, 0, 1] + self.subject_auc[:, 1, 0])
        same_structure = 0.5*(same_map + same_structure_dif_map)
        dif_structure = 0.25*(self.subject_auc[:, 2, 1] + self.subject_auc[:, 2, 0] +\
                              self.subject_auc[:, 2, 1] + self.subject_auc[:, 2, 0])
        main_effect = np.squeeze(same_structure - dif_structure)
        return main_effect



print(os.getcwd())
load_cells = scipy.io.loadmat(os.path.join(os.getcwd(), DATA_PATH, 'AUC_visualisation_EC.mat'))
voxels_data = load_cells['betasAllSubj']
print(voxels_data.shape)
voxels_data = np.transpose(voxels_data, (3, 2, 1, 0))
print(voxels_data.shape)
voxels_ana = Svd_AUC_voxels(voxels_data=voxels_data)
voxels_ana.cal_loo_average()
voxels_ana.run_subjects_auc_calculation()
main_effect = voxels_ana.cal_main_effect()
print(main_effect.shape)
print(np.mean(main_effect, axis=0))

res = stats.ttest_1samp(main_effect, popmean=0, axis=0, nan_policy='propagate', alternative='greater')
print(res.statistic, res.pvalue)

voxels_ana.plot_auc_curve()