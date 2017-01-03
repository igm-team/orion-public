"""
Calculate the Orion score for a given SFS vs. the expected distribution
"""

import numpy as np
import sys
import os
from ConfigParser import ConfigParser
from OrionGlobals import confirm_valid_numerical_argument, get_cfg
from operator import le, lt
from collections import deque

cfg = get_cfg()
effective_population_size_default = cfg.getint(
    "sfs", "effective_population_size")
k_mult_constant = cfg.getfloat("sfs", "k_mult_constant")

class OrionSFS(object):
    """Class representing a __folded__ site frequency spectrum
    Contains data about the sample size, the effective population size,
    the mutation rate in the window, the length of the window, and methods
    for calculating the expected site frequency spectrum under neutrality,
    which are in turn used to calculate the difference from the observed
    in order to calculate the Orion score
    """

    def __init__(self, sample_size, window_length,
                 effective_population_size=effective_population_size_default,
                 ploidy=2, allele_counts=None, sfs=None,
                 mu_vector=None, mu=None):
        """Initialize the SFS; sample_size, effective_population_size, and
        ploidy window_length will all remain constant, whereas the sfs itself,
        and the mutation rate may all change as we consider different regions
        and so do not require specification at initialization
        :param sample_size: the number of individual sequences represented
        :param window_length: the number of bases considered to create the SFS
        :param effective_population_size: the Wright-Fisher effective size of
            the population
        :param ploidy: the number of chromosomes each sequence has
        :param allele_counts: optional deque/list of the allele counts by base
            of the window (insufficiently covered bases
            should be represented by -1)
        :param sfs: optional; the actual site frequency spectrum;
            should be a np.array of length floor(sample_size * ploidy / 2) + 1
            as we include the modification of inclusion
            of the count of non-variant sites
        :param mu_vector: the vector of mus over the window (insufficiently
            covered bases should be represented by -1)
        :param mu: the summed mutation rate over the window
        """
        self.sample_size = confirm_valid_numerical_argument(
            sample_size, "sample_size")
        self.window_length = confirm_valid_numerical_argument(
            window_length, "window_length")
        self.effective_population_size = confirm_valid_numerical_argument(
            effective_population_size, "effective_population_size")
        self.ploidy = confirm_valid_numerical_argument(
            ploidy, "ploidy")
        self.set_xi_and_eta()
        self.num_bases_covered = 0
        if allele_counts:
            if type(allele_counts) is list:
                allele_counts = deque(allele_counts)
            elif type(allele_counts) is not deque:
                raise TypeError("allele_counts must be a deque or list")
            if len(allele_counts) != self.window_length:
                raise ValueError("the length of allele_counts, {l}, does not "
                                 "match the window_length, {w}".format(
                                     l=len(allele_counts),
                                     w=self.window_length))
            if len([x for x in allele_counts if type(x) is not int
                    or x < -1 or x > self.max_eta]):
                raise ValueError("allele_counts contains invalid entries")
            # automatically generate the SFS if the allele counts are provided
            self.sfs = np.zeros(self.max_eta + 1, dtype=np.int64)
            for allele_count in allele_counts:
                if allele_count > -1:
                    self.num_bases_covered += 1
                    self.sfs[allele_count] += 1
            self.allele_counts = allele_counts
            self.set_scoring_weights()
            # the current length of the deque of allele counts
            self.len_window = self.window_length
        elif sfs:
            # just the sfs itself may be specified in order to calculate a score
            # on a static, i.e. non-sliding, window
            self.num_bases_covered = self.window_length
            self.update_sfs(sfs)
            self.set_scoring_weights()
        else:
            self.sfs = np.zeros(self.max_eta + 1, dtype=np.int64)
            self.allele_counts = deque()
            self.len_window = 0
        self.mu = 0.0
        if mu_vector:
            if type(mu_vector) is list:
                mu_vector = deque(mu_vector)
            elif type(mu_vector) is not deque:
                raise TypeError("mu_vector must be a deque or list")
            if len(mu_vector) != self.window_length:
                raise ValueError("the length of mu_vector, {l}, does not "
                                 "match the window_length, {w}".format(
                                     l=len(mu_vector), w=self.window_length))
            for mu_base in mu_vector:
                if mu_base != -1:
                    self.mu += mu_base
            self.mu_vector = mu_vector
        elif mu:
            self.update_mu(mu)
        else:
            self.mu_vector = deque()

        # update fraction of covered bases
        self.update_fraction_of_bases_covered()

    def add_site(self, allele_count, mu_base):
        """Add a site to the window, used for building up the allele_counts
        deque initially
        :param allele_count: the number of minor alleles at this site in the
            sample (specify -1 to indicate site is not covered)
        :param mu_base: the mutation rate for the base
            (specify -1 for site not covered)
        """
        if self.len_window == self.window_length:
            raise ValueError("attempting to add a site to a full length window")
        confirm_valid_numerical_argument(
            allele_count, "allele_count", min_value=-1,
            max_value=self.max_eta, left_op=le)
        self.allele_counts.append(allele_count)
        self.mu_vector.append(mu_base)
        self.len_window += 1
        if allele_count != -1:
            self.num_bases_covered += 1
            self.sfs[allele_count] += 1
            self.mu += mu_base
        if self.len_window == self.window_length:
            # the window is now full sized, so calculate the weights
            self.set_scoring_weights()
            self.update_fraction_of_bases_covered()

    def add_and_pop_site(self, allele_count, mu_base):
        """Add a site to the window on the right, pop a site on the left
        :param allele_count: the number of minor alleles at this site in the
            sample (specify -1 to indicate site is not covered)
        :param mu_base: the mutation rate for the base
            (specify -1 for site not covered)
        """
        if self.len_window != self.window_length:
            raise ValueError("attempting to go to next site when the window "
                             "is not the proper length")
        confirm_valid_numerical_argument(
            allele_count, "allele_count", min_value=-1,
            max_value=self.max_eta, left_op=le)
        old_allele_count = self.allele_counts.popleft()
        old_mu_base = self.mu_vector.popleft()
        self.allele_counts.append(allele_count)
        self.mu_vector.append(mu_base)
        if old_allele_count == -1:
            if allele_count != -1:
                self.increment_bases_covered()
        else:
            self.decrement_eta(old_allele_count)
            self.mu -= old_mu_base
            if allele_count == -1:
                self.decrement_bases_covered()
        if allele_count != -1:
            self.increment_eta(allele_count)
            self.mu += mu_base

    def set_xi_and_eta(self):
        """Initialize max_xi, max_eta, and two vectors for subsequent neutral
        eta estimation that will be continually reused
        """
        # unfolded and folded max allele counts, respectively
        self.max_xi = self.sample_size * self.ploidy
        self.max_eta = self.max_xi / 2
        i = np.arange(1, self.max_eta + 1, dtype=np.float64)
        # k_delta is equal to 0 for all entries, [0, max_alleles - 1] and is
        # equal to 0 if the total possible number of unfolded alleles is odd,
        # otherwise it is 1, so that we do not double count the last element of
        # the SFS
        k_delta = np.zeros(self.max_eta, dtype=np.float64)
        if self.max_xi % 2 == 0:
            k_delta[-1] = 1
        # neutral eta distribution is equal to:
            # theta * (((1 / i) + (1 / (self.max_xi - i))) / (1 + k_delta))
        # so we can just calculate everything other than theta once and store it
        self.neutral_etas = ((1 / i) + (1 / (self.max_xi - i))) / (1 + k_delta)

    def decrement_bases_covered(self):
        """Reduce the number of bases in the window that are covered by 1
        and recalculate the scoring weights
        """
        if self.num_bases_covered == 0:
            raise ValueError("the window has 0 bases covered already and is "
                             "being reduced by 1")
        self.num_bases_covered -= 1
        self.update_fraction_of_bases_covered()

    def increment_bases_covered(self):
        """Increase the number of bases in the window that are covered by 1
        and recalculate the scoring weights
        """
        if self.num_bases_covered == self.window_length:
            raise ValueError("the window already has the maximum number of "
                             "bases covered and is being increased by 1")
        self.num_bases_covered += 1
        self.update_fraction_of_bases_covered()

    def update_fraction_of_bases_covered(self):
        """Update the fraction of bases that are covered in the window
        """
        self.fraction_of_bases_covered = (
            float(self.num_bases_covered) / self.window_length)

    def update_sfs(self, sfs):
        """Update the site frequency spectrum
        :param sfs: the site frequency spectrum as a list or numpy array
        """
        if type(sfs) is list:
            sfs = np.ndarray(sfs, dtype=np.int64)
        elif type(sfs) is not np.ndarray:
            raise TypeError("the SFS must be a list or numpy array")
        if issubclass(sfs.dtype.type, np.integer):
            if sfs.ndim == 1 and sfs.size == (self.max_eta + 1):
                if np.any(sfs < 0):
                    raise ValueError("the SFS has values below 0")
                if sfs.sum() == self.window_length:
                    self.sfs = sfs
                else:
                    raise ValueError("the number of sites in the SFS, {s} "
                                     "does not equal the window length, {l}".
                                     format(s=sfs.sum(), l=self.window_length))
            else:
                raise ValueError("the SFS has invalid dimensions given "
                                 "the sample size and ploidy")
        else:
            raise TypeError("the SFS has a non-integer data type: {dtype}".
                            format(dtype=sfs.dtype.type))

    def update_mu(self, mu):
        """Update the mutation rate for the SFS' window
        :param mu: the mutation rate of the window
        """
        self.mu = confirm_valid_numerical_argument(
            mu, "mu", arg_type=float)

    def increment_eta(self, eta):
        """Increment the count of the number of sites with eta alleles.
        Intended for use in a sliding window where we add one site to the SFS.
        :param eta: the number of alleles for the added site.
        """
        eta = confirm_valid_numerical_argument(
            eta, "eta", min_value=-1, max_value=self.max_eta)
        if self.sfs[eta] < self.window_length:
            self.sfs[eta] += 1
        else:
            raise ValueError("eta_{eta} already has a site count equal to "
                             "the window length".format(eta=eta))

    def decrement_eta(self, eta):
        """Decrement the count of the number of sites with eta alleles.
        Intended for use in a sliding window where we add one site to the SFS.
        :param eta: the number of alleles for the removed site.
        """
        eta = confirm_valid_numerical_argument(
            eta, "eta", left_op=le, max_value=self.max_eta)
        if self.sfs[eta] > 0:
            self.sfs[eta] -= 1
        else:
            raise ValueError("eta_{eta} already has a site count of 0".format(
                eta=eta))

    def get_mu(self):
        """Return mu for the current window.  Intended for use when outputting
        raw data.
        """
        return self.mu
    
    def get_etas(self):
        """Return the vector of etas for the current window.  Intended for use
        when outputting raw data.
        """
        return self.sfs
    
    def calculate_orion_score(self, min_covered_threshold=cfg.get(
        "sfs", "sites_covered_fraction_to_calculate")):
        """Calculate the score for this window.  The primary functionality
        should be left in place but the weighting scheme may be freely
        adjusted/reimplemented.
        :param min_covered_threshold: require that at least this fraction of
            bases are covered to calculate the score, otherwise None is returned
        :return: the Orion score for this window (deviation from neutrality)
        """
        min_covered_threshold = confirm_valid_numerical_argument(
            min_covered_threshold, "min_covered_threshold", arg_type=float,
            max_value=1.0)
        if min_covered_threshold <= self.fraction_of_bases_covered:
            self.sfs_pdf = self.sfs / np.float64(self.num_bases_covered)
            self.neutral_sfs_pdf = self.calculate_neutral_sfs_pdf()
            # confirm valid PDFs
            assert np.isclose(self.sfs_pdf.sum(), 1)
            assert np.isclose(self.neutral_sfs_pdf.sum(), 1)
            return self.calculate_score()
        else:
            return None

    def calculate_neutral_sfs_pdf(self):
        """Calculate the PDF of an equivalent SFS under Wright-Fisher neutrality
        assumptions
        :return: numpy array of the PDF of the neutral SFS
        """
        # theta is equal to the expected number of mutations in the population
        # per generation and the number of mutations along one lineage in one
        # unit of coalescent time
        theta = self._calculate_theta()
        scaled_neutral_etas = theta * self.neutral_etas
        # we prepend the number of bases in the window - the sum of the expected
        # distribution of neutral etas to attain what we'd expect to be the
        # number of non-variant sites
        return (np.concatenate([[self.num_bases_covered - scaled_neutral_etas.sum()],
                               scaled_neutral_etas]) / self.num_bases_covered)

    def set_scoring_weights(self):
        """We weight the deviation from neutrality here, i.e. to prioritize
        signal coming from enrichment of rare variants
        """
        p = np.arange(1, self.max_eta + 2, dtype=np.float64) / self.max_xi
        weights = 1.0 / p
        self.weights = weights / weights.sum()

    def calculate_score(self):
        """Here we take the observed and neutral SFS PDFs, along with the
        weights in order to calculate the score.  This method may be adjusted in
        order to change the formulation of the score.
        :return: the Orion score for this window (deviation from neutrality)
        """
        sfs_pdf_weighted = self.sfs_pdf * self.weights
        neutral_sfs_pdf_weighted = self.neutral_sfs_pdf * self.weights

        # we calculate the score as the weighted mean difference
        score = (sfs_pdf_weighted - neutral_sfs_pdf_weighted).mean()
        
        # divide by theta to normalize for different mutation rates
        score = score / self._calculate_theta()
        
        # multiply by constant for convenience
        score = score * k_mult_constant

        return score

    def _calculate_theta(self):
        """Return the value of theta, which is equal to the expected number
        of mutations in the population per generation and the number of
        mutations along one lineage in one unit of coalescent time.
        :return: Theta.
        """
        return (2.0 * self.effective_population_size * self.ploidy * self.mu)
        
