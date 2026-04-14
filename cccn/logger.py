import logging
import os
import fcntl


class LockedFileHandler(logging.FileHandler):
    """Process-safe file handler using advisory file locks on Linux."""

    def emit(self, record):
        if self.stream is None:
            self.stream = self._open()
        try:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_EX)
            super().emit(record)
        finally:
            fcntl.flock(self.stream.fileno(), fcntl.LOCK_UN)

class Logger:
    def __init__(self, name, level=logging.INFO, log_file=None, rank=None):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        self.logger.propagate = False
        self.logger.handlers.clear()

        rank_text = '' if rank is None else f'[rank={rank}]'
        self.formatter = logging.Formatter(
            f'%(asctime)s - [%(name)s][%(levelname)s]{rank_text}: %(message)s'
        )
        self.console_handler = logging.StreamHandler()
        self.console_handler.setLevel(level)
        self.console_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.console_handler)

        if log_file is None:
            log_file = '{}.log'.format(name)
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)

        # Append to one shared file across MPI ranks.
        self.file_handler = LockedFileHandler(log_file, mode='a')
        self.file_handler.setLevel(level)
        self.file_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.file_handler)

    def get_logger(self):
        return self.logger
    
    def destroy_logger(self):
        self.logger.removeHandler(self.console_handler)
        self.logger.removeHandler(self.file_handler)
        self.file_handler.close()
