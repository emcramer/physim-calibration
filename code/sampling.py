from scipy.stats import qmc, norm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyDOE3 import lhs

class LHSampler:
    """
    A class for performing Latin Hypercube Sampling (LHS) from a parameter dictionary.
    """

    def __init__(self, param_dict: dict, seed:int=123):
        """
        Initialize the LHSampler with a parameter dictionary.

        :param param_dict: A dictionary where each key is a parameter name, and the values are lists of parameter values.
        :param seed: Random seed for the sampler.
        """
        self.param_dict = param_dict
        self.params = list(param_dict.keys())
        self.num_params = len(self.params)
        self.sampler = qmc.LatinHypercube(d=6)  # Initialize with a default dimension
        self.seed = seed
        self.samples=None
        self.raw_samples = None
        self.samples_df = None

    def sample(self, num_samples:int, verbose:bool=False, **kwargs):
        """
        Perform Latin Hypercube Sampling from the parameter dictionary using SciPy.

        :param num_samples: The number of samples to generate.
        :param verbose: Print the samples.
        :param kwargs: Keyword arguments to pass to the SciPy sampler.
        :return: A numpy array of shape (num_samples, len(params)) containing the sampled parameter values.
        """
        # Create a list of value ranges for each parameter
        value_ranges = []
        for param in self.params:
            values = np.array(self.param_dict[param])
            value_ranges.append([values.min(), values.max()])  # Get min and max for scaling

        # Perform Latin Hypercube Sampling using SciPy
        sampler = qmc.LatinHypercube(d=self.num_params, seed=self.seed, **kwargs)  # Set dimension correctly
        samples = sampler.random(num_samples)
        self.raw_samples = samples.copy()
        self.samples = samples

        # Scale and shift samples to match parameter value ranges
        for i, param in enumerate(self.params):
            values = np.array(self.param_dict[param])
            min_val, max_val = value_ranges[i]  # Get min and max for scaling
            self.samples[:, i] = samples[:, i] * (max_val - min_val) + min_val

        if verbose: print(self.samples)

    def get_samples_df(self):
        """
        Returns the sampled parameter values as a pandas data frame.

        :return: A pandas data frame of the sampled parameter values. Column headers are parameter names given in the
                 XML Element Tree format.
        """
        df=pd.DataFrame(
            data=self.samples[0:,0:], 
            index=[i for i in range(self.samples.shape[0])],
            columns=self.params
        )
        self.samples_df = df
        return df

    def set_df_run_ids(self, run_ids):
        """
        Adds the run identifier for each sample as a column in the samples data frame.
        """
        self.samples_df['run_identifier'] = run_ids

    def vis_hist(self):
        """
        Visualize the sampled parameter values for each parameter as histograms.
        """
        for i, param in enumerate(self.params):
            plt.figure()
            plt.hist(self.samples[:, i], bins='auto')
            plt.title(f'Histogram of {param}')
            plt.xlabel(f'{param} value')
            plt.ylabel('Frequency')
            plt.show()
    
    def vis_scatter(self, p1:str, p2:str):
        """
        Visualize the samples in 2D scatter plots.

        :param p1: The name of the first parameter for the scatter plot.
        :param p2: The name of the second parameter for the scatter plot.
        """
        i = list(self.param_dict.keys()).index(p1)
        j = list(self.param_dict.keys()).index(p2)
        plt.figure()
        plt.scatter(self.samples[:, i], self.samples[:, j])
        plt.xlabel(f'{p1}')
        plt.ylabel(f'{p2}')
        plt.show()
    
    def get_discrepancy(self):
        """
        Returns the discrepancy of the of the samples to indicate how well the samples cover the search space.
        :return: The discrepancy of the samples.
        """
        discrepancy = qmc.discrepancy(self.raw_samples)
        return discrepancy

def sampling_sweep(param_dict, n_samples:list=list(range(10, 60, 10)), visualize=True):
    """
    Generates samples for different sample sizes and visualizes the discrepancy of the sampling.
    
    :param param_dict: The parameter dictionary.
    :param n_samples: A list of sample sizes to test.
    :param visualize: Visualize the discrepancy as a function of the sample size.
    """
    discrepancies = [None]*len(n_samples)
    for i, n in enumerate(n_samples):
        sampler = LHSampler(param_dict)
        sampler.sample(n, optimization='random-cd')
        discrepancies[i] = sampler.get_discrepancy()
    print(discrepancies)
    if visualize:
        plt.figure()
        plt.scatter(n_samples, discrepancies)
        plt.plot(n_samples, discrepancies)
        plt.xlabel(f'Samples')
        plt.ylabel(f'Discrepancy')
        plt.title('Coverage of Parameter Search Space')
        plt.show()

##### Tests
def testLHS1():
    """
    Quick unit test for the LHS class.
    """
    param_dict = {'parameter1': [0, 10], 'parameter2': [5, 15], 'parameter3':[20,50]}
    lhs = LHSampler(param_dict)
    samples = lhs.sample(60, verbose=True, optimization='random-cd')
    print(f"Sampling discrepancy: {lhs.get_discrepancy()}")
    lhs.vis_hist()
    lhs.vis_scatter('parameter1', 'parameter2')
    lhs.vis_scatter('parameter1', 'parameter3')
    lhs.vis_scatter('parameter2', 'parameter3')

def testLHS2():
    param_dict = {'parameter1': [0, 10], 'parameter2': [5, 15], 'parameter3':[20,50]}
    sampling_sweep(param_dict, n_samples=list(range(10, 100, 10)))

if __name__ == "__main__":
    #testLHS1()
    testLHS2()