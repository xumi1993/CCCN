import logging

class Logger:
    def __init__(self, name, level=logging.INFO):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)
        self.formatter = logging.Formatter('%(asctime)s - [%(name)s][%(levelname)s]: %(message)s')
        self.console_handler = logging.StreamHandler()
        self.console_handler.setLevel(level)
        self.console_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.console_handler)

    def get_logger(self):
        return self.logger
