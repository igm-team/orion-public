import sys
import os
import re
import argparse

# Need correct configparser based on python version
if sys.version_info[0] == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser

from collections import OrderedDict
import gzip
from operator import le, lt
import subprocess
import shlex
from luigi import Target

cfg = ConfigParser()
cfg.read(os.path.join(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))),
    "cfg/orion.cfg"))
chromosomes = OrderedDict([[chromosome, int(length)] for chromosome, length in
                           cfg.items("chromosomes")])
boundaries_regex = re.compile(cfg.get("atav", "boundaries"))
# http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/
# maintain sort order
CHROMOSOME_LENGTHS = OrderedDict(
    [[chrom, int(length)] for chrom, length in cfg.items("chromosomes")])
vcf_columns = OrderedDict((field, x) for x, field in enumerate(
    cfg.get("ref", "vcf_columns").split(",")))
bases = set(cfg.get("ref", "bases").split(","))

def get_cfg():
    """return the Orion config file
    """
    return cfg

def get_vcf_fields_dict(fields):
    return dict(zip(vcf_columns.keys(), fields))

def get_vcf_fields_line(line):
    """return a dict of the values of a VCF line
    """
    return get_vcf_fields_dict(line.strip().split("\t"))

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    """format help description as specified and print parameter defaults
    """
    pass

def file_exists(file_name):
    """return the real path to the file if it exists, otherwise
    raise an argparse ArgumentTypeError
    """
    if os.path.isfile(file_name):
        return os.path.realpath(file_name)
    else:
        raise argparse.ArgumentTypeError(
            "{file_name} does not exist".format(file_name=file_name))

def float_between_zero_and_one(arg):
    """return the float value of the arg if 0 <= float(arg) <= 1,
    else raise a type error
    """
    try:
        value = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "{arg} is not a valid float".format(arg=arg))
    if value >= 0.0 and value <= 1.0:
        return value
    else:
        raise argparse.ArgumentTypeError(
            "{arg} is not in the range [0, 1]".format(arg=arg))

def non_negative_float(arg):
    """return the float value of the arg if float(arg) >= 0,
    else raise a type error
    """
    try:
        value = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "{arg} is not a valid float".format(arg=arg))
    if value >= 0.0:
        return value
    else:
        raise argparse.ArgumentTypeError(
            "{arg} can not be negative".format(arg=arg))

def natural_number(arg):
    """return the int value of the arg if it's a natural number,
    else raise a type error
    """
    try:
        value = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "{arg} is not a valid int".format(arg=arg))
    if value > 0:
        return value
    else:
        raise argparse.ArgumentTypeError(
            "{arg} is not a valid natural number".format(arg=arg))

def odd_natural_number(arg):
    """return the int value of  the arg if it's a natural number and odd, else
    raise a type error
    """
    value = natural_number(arg)
    if value % 2 == 1:
        return value
    else:
        raise argparse.ArgumentTypeError(
            "{arg} is even".format(arg=arg))

class FormatError(Exception):
    """Thrown if the interval format does not match that expected
    """
    pass

def parse_trinucleotide_mutation_rates(matrix, split_mu_by_mutation=False):
    """verify the passed in matrix meets the defined format
    """
    trinucleotides = {}
    with open(matrix) as matrix_fh:
        header = matrix_fh.next().strip().split("\t")
        try:
            trinucleotide_index = header.index("trinucleotide")
            indexes = {}
            for base in bases:
                indexes[base] = header.index(base)
        except ValueError:
            raise ValueError(
                "{matrix} does not properly define the fields in the header".
                format(matrix=matrix))
        for line in matrix_fh:
            fields = line.strip().split("\t")
            trinucleotide = fields[trinucleotide_index]
            if split_mu_by_mutation:
                for mutation in bases - set(trinucleotide[1]):
                    trinucleotides[(trinucleotide, mutation)] = float(fields[indexes[mutation]])
            else:
                rate = 0.0
                for index in indexes.itervalues():
                    rate += float(fields[index])
                trinucleotides[trinucleotide] = rate
    return trinucleotides

def bed_file(arg):
    """return the parsed BED file regions or raise an exception if it's invalid
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError("{} does not exist".format(arg))

    intervals = [] 
    with open(arg) as bed:
        for line in bed:
            chromosome, start, end = line.strip().split("\t")
            if chromosome in CHROMOSOME_LENGTHS:
                length = CHROMOSOME_LENGTHS[chromosome]
            else:
                raise argparse.ArgumentTypeError(
                    "chromosome {} is not valid".format(chromosome))
            start = int(start)
            end = int(end)
            if 0 <= start < end <= length:
                intervals.append((chromosome, start, end))
            else:
                raise argparse.ArgumentTypeError(
                    "invalid coordinates: {}:{}-{}".format(
                        chromosome, start, end))
    return intervals

def optionally_gzipped_file(arg):
    """return a file handle to the specified argument whether it's gzipped or
    otherwise, and the path to the directory
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(arg + " does not exist")
    dir_name = os.path.dirname(os.path.realpath(arg))
    with open(arg, "rb") as fh:
        magic_number = fh.read(2)
    open_func = gzip.open if magic_number == "\x1f\x8b" else open
    return open_func(arg), dir_name

def directory_exists(arg):
    """return the path to the arg if it's a directory, raise exception otherwise
    """
    rp = os.path.realpath(arg)
    if os.path.isdir(rp):
        return rp
    else:
        raise argparse.ArgumentTypeError("{} does not exist".format(
            rp))

def directory(arg):
    """create the specified directory if necessary and return the real path
    """
    rp = os.path.realpath(arg)
    if not os.path.isdir(rp):
        try:
            os.makedirs(rp)
        except OSError:
            raise argparse.ArgumentTypeError(
                "couldn't create directory {}".format(rp))
    return rp

def create_INFO_dict(INFO):
    """return a dict of values from the INFO field of a VCF
    """
    d = {}
    for entry in INFO.split(";"):
        if "=" in entry:
            variable, value = entry.split("=", 1)
            d[variable] = value
    return d

def mutation_set(arg):
    """return the set of mutations in arg
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError("{} does not exist".format(arg))
    variants = set([])
    with open(arg) as variants_fh:
        for line in variants_fh:
            variants.add(tuple(line.strip().split("\t")))
    return variants

def simplify_REF_ALT_alleles(REF, ALT):
    """take a potentially complex representation of a pair of alleles introduces
    into multiallelic sites and simplify to canonical form
    http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
    """
    # strip shared suffix
    strip = 0
    for x in xrange(1, min(len(REF), len(ALT))):
        if REF[-x] == ALT[-x]:
            strip -= 1
        else:
            break
    if strip:
        REF = REF[:strip]
        ALT = ALT[:strip]
    # strip shared prefix
    strip = 0
    for x in xrange(0, min(len(REF), len(ALT)) - 1):
        if REF[x] == ALT[x]:
            strip += 1
        else:
            break
    # return simplified REF, ALT, and position offset
    return REF[strip:], ALT[strip:], strip

def get_fh(file_name, mode="r"):
    """return a file handle to the file_name whether it's gzipped or not
    """
    dirname = os.path.dirname(file_name)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    if mode == "r":
        with open(file_name) as fh:
            magic_number = fh.read(2)
        if magic_number == "\x1f\x8b":
            return gzip.open(file_name)
        else:
            return open(file_name)
    elif mode == "a":
        if os.path.isfile(file_name):
            with open(file_name) as fh:
                magic_number = fh.read(2)
            if magic_number == "\x1f\x8b":
                return gzip.open(file_name, "a")
        return open(file_name, "a")
    elif mode == "w":
        return (gzip.open(file_name, "w") if file_name.endswith(".gz")
                else open(file_name, "w"))
    else:
        raise ValueError("unknown mode: {}".format(mode))

def mutation_rate_matrix(arg):
    """verify the the passed in matrix meets the defined format
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(
            "{arg} does not exist".format(arg=arg))
    trinucleotides = {}
    with open(arg) as matrix:
        header = matrix.next().strip().split("\t")
        try:
            trinucleotide_index = header.index("trinucleotide")
            A_index = header.index("A")
            C_index = header.index("C")
            G_index = header.index("G")
            T_index = header.index("T")
        except ValueError:
            raise argparse.ArgumentTypeError(
                "{arg} does not properly define the fields in the header".
                format(arg=arg))
        for line in matrix:
            fields = line.strip().split("\t")
            trinucleotide = fields[trinucleotide_index]
            rate = 0.0
            for index in (A_index, C_index, G_index, T_index):
                rate += float(fields[index])
            trinucleotides[codon] = rate
    return trinucleotides

def atav_interval(arg):
    """parse the ATAV formatted interval into genomic coordinates
    """
    m = boundaries_regex.match(arg)
    if not m:
        raise argparse.ArgumentTypeError(
            "{arg} is not in the proper ATAV-defined format".format(arg=arg))
    d = m.groupdict()
    intervals = []
    for interval_pair in d["intervals"].split(","):
        intervals.append(
            [int(position) for position in interval_pair.split("..")])
    return (arg, d["chromosome"], intervals)

def atav_intervals_file(arg):
    """return the parsed ATAV-format intervals in the specified file or raise an
    exception with all that aren't properly formatted
    """
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(
            "{arg} does not exist".format(arg=arg))
    with open(arg) as input_intervals:
        interval_lines = input_intervals.read().splitlines()
    parsed_intervals = []
    misformed_lines = []
    for interval_line in interval_lines:
        m = boundaries_regex.match(interval_line)
        if not m:
            misformed_lines.append(interval_line)
        if misformed_lines:
            continue
        d = m.groupdict()
        intervals = []
        for interval_pair in d["intervals"].split(","):
            intervals.append(
                [int(position) for position in interval_pair.split("..")])
        parsed_intervals.append((interval_line, d["chromosome"], intervals))
    if misformed_lines:
        raise argparse.ArgumentTypeError(
            "the following lines do not match the ATAV-defined format:\n" +
            "".join(misformed_lines))
    return parsed_intervals

def confirm_valid_numerical_argument(
    arg, arg_name, arg_type=int, min_value=0, max_value=sys.maxsize,
    left_op=lt, right_op=le):
    """Confirm that the specified value is valid in the range
        (minimum_value, maximum_value] (by default)
    :param arg: the value to be tested
    :param arg_name: the name of the parameter
    :param arg_type: the type of the parameter, e.g. int or float
    :param min_value: the minimum value for the parameter, exclusive
    :param max_value: the maximum value for the parameter, inclusive
    :param left_op: the operator for testing left_op(min_value, value)
    :param right_op: the operator testing right_op(value, max_value)
    :return: arg_type(arg) if arg is valid
    """
    try:
        value = arg_type(arg)
        if left_op(min_value, value) and right_op(value, max_value):
            return value
        else:
            raise ValueError(
                "{arg_name} ({arg}) is not in the range "
                "{left_endpoint}{min_value}, {max_value}{right_endpoint}".format(
                    arg_name=arg_name, arg=arg, min_value=min_value,
                    max_value=max_value, left_endpoint="(" if left_op == lt else
                    "[", right_endpoint="]" if right_op == le else ")"))
    except TypeError:
        raise TypeError(
            "{arg_name} ({arg}) is not a valid {arg_type}".format(
                arg_name=arg_name, arg=arg, arg_type=arg_type.__name__))

def run_command(command, directory, task_name, print_command=True):
    """Run the specified command in a subprocess and log the output
    """
    if print_command:
        print(command)
    base_name = os.path.join(directory, task_name)
    with open(base_name + ".out", "w") as out_fh, \
            open(base_name + ".err", "w") as err_fh:
        out_fh.write(command + "\n")
        out_fh.flush()
        p = subprocess.Popen(shlex.split(command), stdout=out_fh, stderr=err_fh)
        p.wait()
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, command)

class MultipleFilesTarget(Target):
    """Target for just testing existence of files specified
    """
    def __init__(self, files):
        self.files = files

    def exists(self):
        return all(os.path.isfile(fn) for fn in self.files)

def strip_prefix(string, prefix):
    """Strip the given prefix from the string if it's present
    """
    return string[len(prefix):] if string.startswith(prefix) else string

def strip_suffix(string, suffix):
    """Strip the given suffix from the string if it's present
    """
    return string[:-len(suffix)] if string.endswith(suffix) else string
