import h5py


class NoiseFile:
    def __init__(self, filename, action='w'):
        self.filename = filename
        self.file = h5py.File(self.filename, action)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def write_header(self, header):
        for key, value in header.items():
            self.file.attrs[key] = value
    
    def write_data(self, data, name):
        if name in self.file:
            del self.file[name]
        self.file.create_dataset(name, data=data)