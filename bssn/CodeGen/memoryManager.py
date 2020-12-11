"""
@author: Milinda Fernando(milinda@cs.utah.edu)
@brief: Simple memory manger class for Cache simulation. Only managers the shared main memory, performs simple address translations.
"""


class MemManager():

    def __init__(self,base_addr=0):
        self.__mmap__= dict()
        self.__next_addr = base_addr
    
    def allocate(self,var,num_bytes):
        if var in self.__mmap__:
            print("mem  manager: warning duplicate allocation")

        self.__mmap__[var]=(self.__next_addr,num_bytes)
        self.__next_addr =self.__next_addr + num_bytes
        #print("next add: %d "%self.__next_addr)

    def deallocate(self,var):
        if not (var in self.__mmap__):
            print("mem manager: warning deallocation key does not exit")
            return
        
        del self.__mmap__[var]
    
    def get_address(self,var,k,size_t):
        
        if not (var in self.__mmap__):
            print("mem manager: warning requested key %s does not exist" %var)
            return -1
        
        add_begin = self.__mmap__[var][0]
        add_size  = self.__mmap__[var][1]

        req_add = add_begin + k*size_t
        if(req_add > add_begin + add_size):
            print("mem manager: access out of bounds")

        return req_add




        
    








