class InputFileParameter(Parameter):
    """
    Paramater whose value is an existing file.
    """

    def parse(self, s):
        """
        Returns the real path to the file from the string if it exists.
        """
        if os.path.isfile(s):
            return os.path.realpath(s)
        else:
            raise OSError("{s} does not exist".format(s=s))


class OutputFileParameter(Parameter):
    """
    Parameter whose value is an output file to be created.
    """

    def parse(self, s):
        """
        Returns the real path to the file and attempts to create any
        directories needed.
        """
        s = os.path.realpath(s)
        if not os.path.isdir(os.path.dirname(s)):
            try:
                os.makedirs(os.path.dirname(s))
            except OSError:
                raise OSError(
                  "could not make parent directories for {s}".format(s=s))
        return s
