import logging
from logging.handlers import RotatingFileHandler

class Logger:
    def __init__(self, name, level=logging.INFO):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        self.logger.propagate = False
        self.logger.handlers.clear()
        self.formatter = logging.Formatter('%(asctime)s - [%(name)s][%(levelname)s]: %(message)s')
        self.console_handler = logging.StreamHandler()
        self.console_handler.setLevel(level)
        self.console_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.console_handler)
        # Keep log files bounded in size and preserve only a small history.
        self.file_handler = RotatingFileHandler(
            '{}.log'.format(name),
            mode='w',
            maxBytes=1_000_000,
            backupCount=2,
        )
        self.file_handler.setLevel(level)
        self.file_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.file_handler)

    def get_logger(self):
        return self.logger
    
    def destroy_logger(self):
        self.logger.removeHandler(self.console_handler)
        self.logger.removeHandler(self.file_handler)
        self.file_handler.close()
