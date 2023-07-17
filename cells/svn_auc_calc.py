import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
import os.path
from scipy import interpolate
from scipy.ndimage import gaussian_filter

class Svn_AUC:
    def __init__(self, grid_cells, place_cells):
        self.grid_cells = grid_cells
        self.place_cells = place_cells