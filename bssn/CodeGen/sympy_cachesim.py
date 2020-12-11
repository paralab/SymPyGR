"""
@author: Milinda Fernando (milinda@cs.utah.edu)
School of Computing, University of Utah.
@brief: Cache simulator for evaluating sympy expressions
"""

import cachesim


"""
Simulate simple cache hierarchy to evaluate sympy expressions
"""
class SympyCacheSim:
    
    def __init__(self):
        self.__mem__ = cachesim.MainMemory()
        self.__l3__  = cachesim.Cache("L3", 20480, 16, 64, "LRU")  # 20MB: 20480 sets, 16-ways with cacheline size of 64 bytes
        
        self.__mem__.load_to(self.__l3__)
        self.__mem__.store_from(self.__l3__)
        
        self.__l2__ = cachesim.Cache("L2", 256, 8, 64, "LRU", store_to=self.__l3__, load_from=self.__l3__)  # 256KB
        self.__l1__ = cachesim.Cache("L1", 32, 8, 64, "LRU", store_to=self.__l2__, load_from=self.__l2__)  # 32KB
        self.__cs__ = cachesim.CacheSimulator(self.__l1__, self.__mem__)
        
    def get_cachesim(self):
        return self.__cs__


    
    




