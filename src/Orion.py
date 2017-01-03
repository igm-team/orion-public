"""
Calculate Orion scores genome-wide or in intervals from a BED file
with cluster parallelization using luigi
Windows are (arbitrarily) processed in windows of size 100,000
"""

import luigi
from luigi.contrib.sge import SGEJobTask
import os
from ConfigParser import ConfigParser
from OrionGlobals import *
from OrionSFS import OrionSFS
from operator import le, lt
import tabix
from collections import OrderedDict
import subprocess
from sys import maxsize
from MD5Target import MD5Target

cfg = get_cfg()
chromosomes = OrderedDict([[chromosome, int(length)] for chromosome, length in
                           cfg.items("chromosomes")])

def parse_intervals(bed_file):
    """Return the list of intervals specified in the BED file
    :param bed_file: BED file of intervals to parse
    :return: list of intervals from the BED file
    """
    if not os.path.isfile(bed_file):
        raise OSError("{bed_file} does not exist".format(bed_file=bed_file))
    intervals = []
    with open(bed_file) as bed:
        for line in bed:
            line = line.strip()
            chromosome, start, end = line.split("\t")
            start = int(start)
            end = int(end)
            if 0 <= start < end < chromosomes[chromosome]:
                intervals.append((chromosome, start, end))
            else:
                raise ValueError(
                    "interval [{start}, {end}) outside of chromosome "
                    "{chromosome} range [0, {c_end})".format(
                        start=start, end=end, chromosome=chromosome,
                        c_end=chromosomes[chromosome]))
    return intervals


class OrionTask(SGEJobTask):
    """Class for calculating Orion scores in a window on the cluster
    :param sample-size: the number of sequences in the sample
    :param ploidy: the ploidy of the sequences
    :param effective-population-size: the effective population size
    :param coverage-file: coverage file per position
    :param allele-counts-file: allele counts per variant site file
    :param genome-mu-file: genome-wide mutation rate file
    :param samples-covered-fraction: the fraction of samples that must have at
        least one read for a site to be considered
    :param mean-coverage: the minimum average coverage for samples with
        some coverage for a site to be considered
    :param minimum-sites-covered-fraction: the minimum fraction of bases that
        must be covered at the minimum thresholds in order to calculate a score
    :param window-length: the length of the SFS window
    :param chromosome: the chromosome to calculate scores for
    :param start: the position to start calculating scores
    :param end: the position to end calculating scores
    :param output-directory: the directory to output to
    :param full-check: test MD5 checksums of intermediate files instead of just
         testing files' existence
    :param raw-data: output the raw data instead of calculating scores
    """
    sample_size = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=cfg.getint("sfs", "sample_size"),
        description="the number of sequences in the sample")
    ploidy = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=2, description="the ploidy of the sequences")
    effective_population_size = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int,
        default=cfg.getint("sfs", "effective_population_size"),
        description="the effective size of the population")
    coverage_file = luigi.InputFileParameter(
        default=cfg.get("sfs", "coverage"),
        description="coverage file per position")
    allele_counts_file = luigi.InputFileParameter(
        default=cfg.get("sfs", "allele_counts"),
        description="allele counts per variant site file")
    genome_mu_file = luigi.InputFileParameter(
        default=cfg.get("ref", "genome_mu"),
        description="the genome-wide mutation rate file")
    samples_covered_fraction = luigi.NumericalParameter(
        min_value=0, max_value=1, var_type=float, left_op=lt, right_op=le,
        default=cfg.getfloat("global", "min_samples_cvg_fraction"),
        description="the minimum fraction of samples that have to have "
        "coverage for a site to be considered")
    minimum_sites_covered_fraction = luigi.NumericalParameter(
        min_value=0, max_value=1, var_type=float, left_op=lt, right_op=le,
        default=cfg.getfloat("sfs", "sites_covered_fraction_to_calculate"),
        description="the minimum fraction of bases that must be covered at the"
        " minimum thresholds in order to calculate a score")
    window_length = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=cfg.getint("sfs", "window_length"),
        description="the length of the SFS window")
    chromosome = luigi.ChoiceParameter(
        choices=chromosomes.keys(),
        description="the chromosome to calculate scores for")
    start = luigi.NumericalParameter(
        min_value=0, max_value=maxsize, var_type=int,
        description="the position to start calculating scores")
    end = luigi.NumericalParameter(
        min_value=2, max_value=maxsize, var_type=int,
        description="the position to end calculating scores") 
    output_directory = luigi.Parameter(
        description="the directory to output to")
    shared_tmp_dir = "/nfs/seqscratch09/orion_jobs"
    full_check = luigi.BoolParameter(
        significant=False,
        description="test MD5 checksums of intermediate files instead of "
        "just testing files' existence (adds a fair bit of overhead)")
    raw_data = luigi.BoolParameter(
        description="output the raw data instead of calculating scores")

    def __init__(self, *args, **kwargs):
        super(OrionTask, self).__init__(*args, **kwargs)
        self.half_window_length = self.window_length / 2

    def work(self):
        pos_after_decimal_to_round = cfg.getint(
            "sfs", "pos_after_decimal_to_round")
        output_fh = self.output().open("w")
        orion_sfs = OrionSFS(
            sample_size=self.sample_size,
            effective_population_size=self.effective_population_size,
            window_length=self.window_length, ploidy=self.ploidy)
        cov = tabix.open(self.coverage_file)
        var = tabix.open(self.allele_counts_file)
        mu = tabix.open(self.genome_mu_file)
        # indices are inclusive in tabix, whereas we have an exclusive end
        cov_iter = cov.query(self.chromosome, self.start, self.end - 1)
        var_iter = var.query(self.chromosome, self.start, self.end - 1)
        mu_iter = mu.query(self.chromosome, self.start, self.end - 1)
        start_output_pos = self.start + self.window_length - 1
        last_position = self.start - 1
        try:
            next_pos, next_count = [int(value) for value in var_iter.next()[1:]]
            var_done = False
        except StopIteration:
            var_done = True
        for (_, position, samples_covered_fraction,
             samples_covered_above_20gq_fraction, mean_coverage) in cov_iter:
            position = int(position)
            if position != (last_position + 1):
                raise ValueError("incorrect position in coverage file; expected"
                                 "{expected}, got {actual}".format(
                                     expected=last_position + 1,
                                     actual=position))
            last_position = position
            samples_covered_fraction = float(samples_covered_fraction)
            samples_covered_above_20gq_fraction = float(
                samples_covered_above_20gq_fraction)
            mean_coverage = float(mean_coverage)
            _, mu_position, mu_base = mu_iter.next()
            mu_position = int(mu_position)
            mu_base = float(mu_base)
            if mu_position != position:
                raise ValueError("incorrect position in mutation rate file; "
                                 "expected {expected}, got {actual}".format(
                                     expected=position, actual=mu_position))
            if not var_done and next_pos == position:
                allele_count = next_count
                try:
                    next_pos, next_count = [
                        int(value) for value in var_iter.next()[1:]]
                except StopIteration:
                    var_done = True
            else:
                allele_count = 0
            if samples_covered_above_20gq_fraction < self.samples_covered_fraction:
                allele_count = -1
                mu_base = -1
            if position <= start_output_pos:
                orion_sfs.add_site(allele_count, mu_base)
            else:
                orion_sfs.add_and_pop_site(allele_count, mu_base)
            if position >= start_output_pos:
                out_pos = position - self.half_window_length
                if self.raw_data:
                    etas = orion_sfs.get_etas()
                    etas_list = []
                    for x, eta in enumerate(etas):
                        if eta:
                            etas_list.append("{x}:{eta}".format(x=x, eta=eta))
                    etas_string = " ".join(etas_list)
                    output_fh.write(
                        "{chromosome}\t{pos}\t{fraction_covered}\t{mu}\t{etas_string}\n".
                        format(
                            chromosome=self.chromosome, pos=out_pos,
                            fraction_covered=round(
                                orion_sfs.fraction_of_bases_covered, 2),
                            mu=orion_sfs.get_mu(), etas_string=etas_string))
                else:
                    orion_score = orion_sfs.calculate_orion_score(
                        min_covered_threshold=self.minimum_sites_covered_fraction)
                    output_fh.write(
                        "{chromosome}\t{pos}\t{pos_plus_1}\t{score}\t{fraction_covered}\n".
                        format(
                            chromosome=self.chromosome,
                            pos=out_pos, pos_plus_1=out_pos + 1,
                            score="NA" if orion_score is None else
                            round(orion_score, pos_after_decimal_to_round),
                            fraction_covered=round(
                                orion_sfs.fraction_of_bases_covered, 2)))
        output_fh.close()

        # verify the output and MD5 checksum it
        error_msg = "failed verifying output"
        with self.output().open() as in_file:
            for line_num, line in enumerate(in_file):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise ValueError("incorrect number of fields @line #{} in {}".format(
                        line_num, in_file.name))
                if fields[0] != self.chromosome:
                    raise ValueError(error_msg)
                if int(fields[1]) != (
                    line_num + self.start + self.half_window_length):
                    raise ValueError(
                        error_msg + ": incorrect position @ line #{} in {}".
                        format(line_num, in_file.fn))
        expected_line_num = self.end - self.start - self.window_length
        if line_num != expected_line_num:
            # invalid number of lines
            raise ValueError(
                error_msg + "; lines:{line_num}; expected:{expected}".format(
                    line_num=line_num, expected=expected_line_num))
        cmd = "md5sum {out_file} > {out_file}.md5".format(
            out_file=self.output().outputs()[0] if self.full_check else
            self.output().fn)
        subprocess.check_call(cmd, shell=True)

    def output(self):
        out_file = os.path.join(
            self.output_directory, "{chromosome}.{start}.{end}.txt".format(
                chromosome=self.chromosome, start=self.start, end=self.end))
        if self.full_check:
            return MD5Target(files=[out_file])
        else:
            return luigi.LocalTarget(out_file)

class CalculateOrionScores(luigi.Task):
    """Primary class for calculating Orion scores in a sliding window
    :param bed-file: optional BED file to limit to regions of interest
    :param sample-size: the number of sequences in the sample
    :param ploidy: the ploidy of the sequences
    :param effective-population-size: the effective population size
    :param coverage-file: coverage file per position
    :param allele-counts-file: allele counts per variant site file
    :param genome-mu-file: genome-wide mutation rate file
    :param samples-covered-fraction: the fraction of samples that must have at
        least one read for a site to be considered
    :param mean-coverage: the minimum average coverage for samples with
        some coverage for a site to be considered
    :param minimum-sites-covered-fraction: the minimum fraction of bases that
        must be covered at the minimum thresholds in order to calculate a score
    :param window-length: the length of the SFS window
    :param bases-per-job: the number of bases per task to run
    :param run-locally: process the batches of scores locally instead of on the
        cluster
    :param output_path: the output file location
    :param output-directory: the output directory for temporary files
    :param full-check: test MD5 checksums of intermediate files instead of just
         testing files' existence
    :param raw-data: output the raw data instead of calculating scores
    """
    bed_file = luigi.InputFileParameter(
        default=None,
        description="an optional BED file to restrict score calculation to")
    sample_size = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=cfg.getint("sfs", "sample_size"),
        description="the number of sequences in the sample")
    ploidy = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=2, description="the ploidy of the sequences")
    effective_population_size = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int,
        default=cfg.getint("sfs", "effective_population_size"),
        description="the effective size of the population")
    coverage_file = luigi.InputFileParameter(
        default=cfg.get("sfs", "coverage"),
        description="coverage file per position")
    allele_counts_file = luigi.InputFileParameter(
        default=cfg.get("sfs", "allele_counts"),
        description="allele counts per variant site file")
    genome_mu_file = luigi.InputFileParameter(
        default=cfg.get("ref", "genome_mu"),
        description="the genome-wide mutation rate file")
    samples_covered_fraction = luigi.NumericalParameter(
        min_value=0, max_value=1, var_type=float, left_op=lt, right_op=le,
        default=cfg.getfloat("global", "min_samples_cvg_fraction"),
        description="the minimum fraction of samples that have to have "
        "coverage for a site to be considered")
    minimum_sites_covered_fraction = luigi.NumericalParameter(
        min_value=0, max_value=1, var_type=float, left_op=lt, right_op=le,
        default=cfg.getfloat("sfs", "sites_covered_fraction_to_calculate"),
        description="the minimum fraction of bases that must be covered at the"
        " minimum thresholds in order to calculate a score")
    window_length = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=cfg.getint("sfs", "window_length"),
        description="the length of the SFS window")
    bases_per_job = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=cfg.getint("sfs", "bases_per_job"), significant=False,
        description="the number of bases per task to run")
    run_locally = luigi.BoolParameter(
        description="process the batches of scores locally", significant=False)
    poll_time = luigi.NumericalParameter(
        min_value=1, max_value=maxsize, var_type=int, right_op=le,
        default=10, description="the time to wait before checking the status "
        "of SGE cluster jobs", significant=False)
    dont_remove_tmp_dir = luigi.BoolParameter(
        significant=False, description="don't delete the temporary directories "
        "of the files submitted to the cluster")
    output_path = luigi.OutputFileParameter(description="the output file location")
    output_directory = luigi.Parameter(
        description="the output directory for temporary files")
    full_check = luigi.BoolParameter(
        significant=False,
        description="test MD5 checksums of intermediate files instead of "
        "just testing files' existence (adds a fair bit of overhead)")
    raw_data = luigi.BoolParameter(
        description="output the raw data instead of calculating scores")

    def __init__(self, *args, **kwargs):
        super(CalculateOrionScores, self).__init__(*args, **kwargs)
        if self.bed_file:
            raw_intervals = parse_intervals(self.bed_file)
        else:
            raw_intervals = [(
                chromosome, 0, end_position) for chromosome, end_position
                in chromosomes.iteritems()]
        self.half_window_length = self.window_length / 2
        if not self.output_path.endswith(".gz"):
            self.output_path = self.output_path + ".gz"
        self.intervals = []
        # create intervals with overlap of window_length / 2 for each raw
        # interval in order to process in parallel
        for chromosome, start, end in raw_intervals:
            adjusted_start = max(0, start - self.half_window_length)
            adjusted_end = min(
                chromosomes[chromosome], end + self.half_window_length)
            for new_start in xrange(
                adjusted_start, adjusted_end, self.bases_per_job + 1):
                self.intervals.append((
                    chromosome, max(0, new_start),
                    min(new_start + self.window_length + self.bases_per_job,
                        adjusted_end, chromosomes[chromosome])))

    def requires(self):
        return [self.clone(
            OrionTask, chromosome=chromosome, start=start, end=end)
            for chromosome, start, end in self.intervals]

    def run(self):
        """Concatenate the results of each OrionTask run on a subset of the data
        """
        header = (cfg.get("sfs", "orion_header_raw") if self.raw_data
                      else cfg.get("sfs", "orion_header"))
        header = header.replace(",", "\t")
        with open(os.path.splitext(self.output_path)[0], "w") as out:
            out.write(header + "\n")
            for input_target in self.input():
                in_fh = (open(input_target.outputs()[0]) if self.full_check
                         else open(input_target.fn))
                for line in in_fh:
                    out.write(line)
                in_fh.close()

        with open(os.path.splitext(self.output_path)[0]) as in_file:
            header_in_file = in_file.next().strip()
            if header != header_in_file:
                return False
            error_msg = "failed verifying output"
            line_count = 0
            for chromosome, start, end in self.intervals:
                for position in xrange(
                    start + self.half_window_length, end -
                    self.half_window_length):
                    try:
                        next_line = in_file.next().rstrip("\n")
                        fields = next_line.split("\t")
                        if len(fields) != 5:
                            raise ValueError(error_msg)
                        chrom, pos, = fields[:2]
                        if chromosome != chrom:
                            raise ValueError(error_msg)
                        if int(pos) != position:
                            raise ValueError(error_msg)
                        if self.raw_data:
                            covered, mu, etas = fields[2:]
                            # just verify appropriate types here
                            float(mu)
                            [[int(value) for value in eta.split(":")]
                             for eta in etas.split(" ") if eta]
                        else:
                            pos_plus_one, score, covered = fields[2:]
                            if (int(pos) + 1) != int(pos_plus_one):
                                raise ValueError(
                                    error_msg + " at line {}".format(line_count))
                            if score != "NA":
                                float(score)
                        float(covered)
                    except:
                        raise ValueError(
                            error_msg + " at line {}".format(line_count))
                    line_count += 1
            try:
                in_file.next()
                # if this worked, we didn't reach the end of the file
                raise ValueError(error_msg + ": file too long")
            except StopIteration:
                pass
        subprocess.check_call(["bgzip", os.path.splitext(self.output_path)[0]])
        subprocess.check_call(["tabix", "-s", "1", "-b", "2", "-e", "2",
                               self.output_path])
        subprocess.check_call("md5sum {orion} > {orion}.md5".format(
            orion=self.output_path), shell=True)
        subprocess.check_call("md5sum {orion_tbi} > {orion_tbi}.md5".format(
            orion_tbi=self.output_path + ".tbi"), shell=True)
    
    def output(self):
        return MD5Target(files=[self.output_path, self.output_path + ".tbi"])

if __name__ == "__main__":
    luigi.run()
