import luigi
import subprocess
import os
from OrionGlobals import get_fh

class MD5Target(luigi.Target):
    """A Target that verifies that the MD5 checkum(s) of the file(s) match
    the .md5 file(s)
    :param files: a list of files to check
    :param full_check: check the MD5 checksum(s) (otherwise just check file
        existence)
    """
    def __init__(self, files, full_check=True):
        self.files = files
        self.full_check = full_check

    def exists(self):
        for target_file in self.files:
            if not os.path.isfile(target_file):
                return False
            if not os.path.isfile(target_file + ".md5"):
                return False
        if self.full_check:
            for target_file in self.files:
                # only re-run if the MD5 checksum doesn't match
                # the previous results
                if subprocess.call(["md5sum", "-c", target_file + ".md5"]):
                    return False
        return True

    def open(self, mode="r"):
        return get_fh(self.files[0], mode)

    def outputs(self):
        return self.files
