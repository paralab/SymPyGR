"""
Cache simulation for the BSSNKO equations
"""

import bssnDerivs
import bssn_test as bssn
import memoryManager 
import sympy_cachesim 
import dendroutils


"""
construct a memory map for single block.
 blk_size : The block is 3d and the blk_size=[nx,ny,nz] sizes of the each dimension
 mm : memory manager
"""
def single_blk_malloc(blk_size,mm,size_t,suffix="[pp]"):
    
    n=(blk_size[0]*blk_size[1]*blk_size[2])

    # BSSNKO in variables. 
    for d in bssnDerivs.BSSN_IN:
        mm.allocate(d+suffix,size_t*n)

    for d in bssnDerivs.BSSN_OUT:
        mm.allocate(d+suffix,size_t*n)
    
    # first derivs
    for d in bssnDerivs.FUNC_D_I:
        mm.allocate(d+suffix,size_t*n)
    
    # second derivs
    for d in bssnDerivs.FUNC_D_IJ:
        mm.allocate(d+suffix,size_t*n)

    # add the adv derivs (we are not currently using this)
    for d in bssnDerivs.FUNC_AD_I:
        mm.allocate(d+suffix,size_t*n)
    
    return



"""
driver 1 : 
for each expr: 
    for each pt in blk:
        do computation. 
"""

def cachesim_driver1(expr_dict, mm, cache, size_t, iter_space, blk_sz):
    
    
    print("driver 1 : for each expr : for each grid_pt")

    dep_set_pp = dict()
    for (v,e) in expr_dict.items():
        dep_set = dendroutils.advanced_free_symbols(e)

        dep_set_pp_list=list()
        for d in dep_set:
            if "[pp]" in str(d):
                dep_set_pp_list.append(d)

        dep_set_pp[v]=dep_set_pp_list
        #print("expr %s : \n  dep set_pp : %s" %(v,dep_set_pp[v]))

        
        
    for (v,e) in expr_dict.items():
        for k in iter_space[2]:
            for j in iter_space[1]:
                for i in iter_space[0]:
                    pp = k*blk_sz[0]*blk_sz[1] + j*blk_sz[0] + i
                    address_rhs = mm.get_address(str(v),pp,size_t)
                    cache.load(address_rhs,length=size_t)
                    for d in dep_set_pp[v]:
                        dep_address = mm.get_address(str(v),pp,size_t)
                        cache.load(dep_address,length=size_t)

                    cache.store(address_rhs,length=size_t)

def cachesim_driver2(expr_dict, mm, cache, size_t, iter_space, blk_sz):
    
    print("driver 2 : for each grid_pt : for each expr")
    
    dep_set_pp = dict()
    dep_unique=set()
    for (v,e) in expr_dict.items():
        dep_set = dendroutils.advanced_free_symbols(e)

        dep_set_pp_list=list()
        for d in dep_set:
            if "[pp]" in str(d):
                dep_set_pp_list.append(d)

        
        dep_unique.update(dep_set_pp_list)
        dep_set_pp[v]=dep_set_pp_list
        #print("expr %s : \n  dep set_pp : %s" %(v,dep_set_pp[v]))
        print("BSSN total [pp] dep: %d" %len(dep_unique))
        
    for k in iter_space[2]:
        for j in iter_space[1]:
            for i in iter_space[0]:
                pp = k*blk_sz[0]*blk_sz[1] + j*blk_sz[0] + i
                for (v,e) in expr_dict.items():
                    address_rhs = mm.get_address(str(v),pp,size_t)
                    cache.load(address_rhs,length=size_t)
                    for d in dep_set_pp[v]:
                        dep_address = mm.get_address(str(v),pp,size_t)
                        cache.load(dep_address,length=size_t)
                    cache.store(address_rhs,length=size_t)                        


################
# parameters

PW=3
blk_lev=0
ele_order=6

blk_sz_1d = (ele_order+1) + ( (1<<blk_lev) -1 )*ele_order + 2*PW
blk_sz =[blk_sz_1d,blk_sz_1d,blk_sz_1d]
iter_space = [range(PW,blk_sz[0]-PW),range(PW,blk_sz[1]-PW),range(PW,blk_sz[2]-PW)]
print("parameters eleOrder : %d blk_sz: %s iter_range: %s " %(ele_order,blk_sz,iter_space))

################

mem_manager = memoryManager.MemManager(0)
single_blk_malloc(blk_sz,mem_manager,8)

#print(mem_manager.__mmap__)
#for (k,v) in mem_manager.__mmap__.items():
#    print("var %s value: %s" %(k,v))

cache       = sympy_cachesim.SympyCacheSim().get_cachesim()
expr_dict   = dendroutils.extract_expressions(bssn.outs,bssn.vnames)

print("warm up run begin")
## cache warm up run
for w_iter in range(0,2):
    cachesim_driver2(expr_dict,mem_manager,cache,8,iter_space,blk_sz)
    cache.force_write_back()

cache.reset_stats()

print("warm up run end")


print("Actual run begin")
cachesim_driver2(expr_dict,mem_manager,cache,8,iter_space,blk_sz)
cache.print_stats()



#cache1 = sympy_cachesim.SympyCacheSim().get_cachesim()
#cache1.load(0,length=1)
#cache1.load(1,length=1)
#cache1.load(8,length=1)
#cache1.load(9,length=1)
#cache1.force_write_back()

#cache1.print_stats()





