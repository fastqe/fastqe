import logging
import sys

def init_logging(loglevel="DEBUG", logfile=None):
    """
    Initializes the logging system with a specified log level and custom formatter.

    Parameters:
    loglevel (str): The log level to use (e.g., "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL").
                    Defaults to "DEBUG".
    logfile (str, optional): The file to which log messages should be written. If None, log messages
                             will only be output to stderr. Defaults to None.

    Raises:
    ValueError: If the provided log level is invalid.
    """
    # Convert log level string to numeric value
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {loglevel}")

    # Define log message format
    fmt = "%(asctime)s | %(levelname)s | %(message)s"

    # Create a StreamHandler to output log messages to stderr
    handler_sh = logging.StreamHandler(sys.stderr)
    handler_sh.setFormatter(logging.Formatter(fmt))

    # Configure the logging system
    if logfile is not None:
        # Create a FileHandler to output log messages to a file
        handler_fh = logging.FileHandler(logfile)
        handler_fh.setFormatter(logging.Formatter(fmt))
        logging.basicConfig(format=fmt, level=numeric_level, handlers=[handler_sh, handler_fh])
    else:
        logging.basicConfig(format=fmt, level=numeric_level, handlers=[handler_sh])

class CustomFormatter(logging.Formatter):
    """
    Custom logging formatter to add color to log messages based on their severity level.
    """

    # Adapted from https://alexandra-zaharia.github.io/posts/make-your-own-custom-color-formatter-with-python-logging
    # ANSI escape codes for colors
    grey = "\x1b[38;21m"
    blue = "\x1b[38;5;39m"
    yellow = "\x1b[38;5;226m"
    red = "\x1b[38;5;196m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    def __init__(self, fmt):
        """
        Initializes the CustomFormatter with a specified format string.

        Parameters:
        fmt (str): The format string for log messages.
        """
        super().__init__()
        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset,
        }

    def format(self, record):
        """
        Formats a log record with the appropriate color based on its severity level.

        Parameters:
        record (logging.LogRecord): The log record to format.

        Returns:
        str: The formatted log message.
        """
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
