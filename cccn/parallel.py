import mpi4py.MPI as MPI
import numpy as np


class MyMPI(object):
    def __init__(self):
        self.world_comm = MPI.COMM_WORLD
        self.world_rank = self.world_comm.Get_rank()
        self.world_size = self.world_comm.Get_size()
        self.node_comm = self.world_comm.Split_type(MPI.COMM_TYPE_SHARED)
        self.node_rank = self.node_comm.Get_rank()
        self.node_size = self.node_comm.Get_size()
        self.node_rank_map = np.zeros(self.world_size, dtype=int)
        self.node_rank_map[self.world_rank] = self.node_rank
        self.node_rank_map = self.sum_all_all(self.node_rank_map)

    def prepare_shm(self, shape:list, dtype):
        itemsize = np.dtype(dtype).itemsize
        n_elem = np.prod(shape)
        n = n_elem if self.node_rank == 0 else 0
        win = MPI.Win.Allocate_shared(n * itemsize, itemsize, comm=self.node_comm)
        
        buf, itemsize = win.Shared_query(0)
        buffer = np.ndarray(buffer=buf, dtype=dtype, shape=shape)
        
        win.Fence()
        return buffer, win
    
    def sync_from_main(self, data):
        if (self.world_rank == 0):
            for i in range(1, self.world_size):
                if self.node_rank_map[i] == 0:
                    self.send(data, dest=i)
        elif self.node_rank == 0:
            data = self.recv(source=0)
    
    def finalize_mpi(self):
        self.world_comm.Barrier()
        self.world_comm.free()
        self.node_comm.free()

    def synchronize_all(self):
        self.world_comm.Barrier()

    def bcast(self, data):
        return self.world_comm.bcast(data, root=0)

    def sum_all(self, data):
        return self.world_comm.reduce(data, op=MPI.SUM, root=0)
    
    def sum_all_all(self, data):
        return self.world_comm.allreduce(data, op=MPI.SUM)

    def max_all(self, data):
        return self.world_comm.reduce(data, op=MPI.MAX, root=0)
    
    def max_all_all(self, data):
        return self.world_comm.allreduce(data, op=MPI.MAX)
    
    def min_all(self, data):
        return self.world_comm.reduce(data, op=MPI.MIN, root=0)
    
    def min_all_all(self, data):
        return self.world_comm.allreduce(data, op=MPI.MIN)
    
    def send(self, data, dest):
        self.world_comm.send(data, dest=dest)

    def recv(self, source):
        return self.world_comm.recv(source=source)
    
    def select_index(self, n):
        idx = []
        for i in range(n):
            dist_rank = i % self.world_size
            if dist_rank == self.world_rank:
                idx.append(i)
        return idx
    



    
        

    
