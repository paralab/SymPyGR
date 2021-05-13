"""
@author: Milinda Fernando (milinda@cs.utah.edu)
School of Computing, University of Utah.
@brief: Cache simulator for evaluating sympy expressions

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use,copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
"""

import cachesim
"""
(base) [xxxxxxx@kingspeak1 ~]$ cat kp_lscpu
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                56
On-line CPU(s) list:   0-55
Thread(s) per core:    2
Core(s) per socket:    14
Socket(s):             2
NUMA node(s):          2
Vendor ID:             GenuineIntel
CPU family:            6
Model:                 79
Model name:            Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
Stepping:              1
CPU MHz:               2317.529
CPU max MHz:           3300.0000
CPU min MHz:           1200.0000
BogoMIPS:              4794.77
Virtualization:        VT-x
L1d cache:             32K
L1i cache:             32K
L2 cache:              256K
L3 cache:              35840K
NUMA node0 CPU(s):     0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,
                       30,32,34,36,38,40,42,44,46,48,50,52,54
NUMA node1 CPU(s):     1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,
                       31,33,35,37,39,41,43,45,47,49,51,53,55
Flags:                 fpu vme de pse tsc msr pae mce cx8 apic sep mtrr
                       pge mca cmov pat pse36 clflush dts acpi mmx fxsr
                       sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm
                       constant_tsc arch_perfmon pebs bts rep_good nopl
                       xtopology nonstop_tsc aperfmperf eagerfpu pni
                       pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2
                       ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1
                       sse4_2 x2apic movbe popcnt aes xsave avx f16c
                       rdrand lahf_lm abm 3dnowprefetch epb cat_l3 cdp_l3
                       intel_pt tpr_shadow vnmi flexpriority ept vpid fsgsbase
                       tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm
                       rdt_a rdseed adx smap xsaveopt cqm_llc cqm_occup_llc
                       cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts
"""
"""
Simulate simple cache hierarchy to evaluate sympy expressions
"""


class SympyCacheSim:
    def __init__(self):
        self.__mem__ = cachesim.MainMemory()
        # 20MB: 20480 sets, 16-ways with cacheline size of 64 bytes
        self.__l3__ = cachesim.Cache("L3", 20480, 16, 64, "LRU")

        self.__mem__.load_to(self.__l3__)
        self.__mem__.store_from(self.__l3__)

        # NOTE: these values aren't the same as the example
        self.__l2__ = cachesim.Cache("L2",
                                     256,
                                     8,
                                     64,
                                     "LRU",
                                     store_to=self.__l3__,
                                     load_from=self.__l3__)  # 256KB
        self.__l1__ = cachesim.Cache("L1",
                                     32,
                                     8,
                                     64,
                                     "LRU",
                                     store_to=self.__l2__,
                                     load_from=self.__l2__)  # 32KB
        self.__cs__ = cachesim.CacheSimulator(self.__l1__, self.__mem__)

    def get_cachesim(self):
        return self.__cs__


if __name__ == "__main__":

    # quick test of the cache for reference
    # this was taken right from PyCacheSim's documentation
    test_cache = SympyCacheSim()

    # laod the cache as set from before
    cs = test_cache.get_cachesim()

    # load a byte from address 2342, should be a miss in all cache
    cs.load(2342)
    # store 8 bytes to addresses 512-519, should also be a load miss
    cs.store(512, length=8)
    # load address 512 until (exclusive) 520 (eight bytes)
    cs.load(512, length=8)

    # force write back, then print stats
    cs.force_write_back()
    cs.print_stats()

    # NOTE: in the print statement:
    # each row referrs to a single memory level, starting with L1 and ending
    # with main Mem
    # Then we have the information about how many hits were made in each type,
    # then misses, then loads, then store, and then evictions
    # take a look at their official GitHub for more info
    # https://github.com/RRZE-HPC/pycachesim
