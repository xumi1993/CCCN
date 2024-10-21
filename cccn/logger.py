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
        self.file_handler = logging.FileHandler('{}.log'.format(name))
        self.file_handler.setLevel(level)
        self.file_handler.setFormatter(self.formatter)
        self.logger.addHandler(self.file_handler)

    def get_logger(self):
        return self.logger
    
    def destroy_logger(self):
        self.logger.removeHandler(self.console_handler)
