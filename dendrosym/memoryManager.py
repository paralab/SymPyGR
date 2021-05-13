"""
@author: Milinda Fernando(milinda@cs.utah.edu)
@brief: Simple memory manger class for Cache simulation. Only managers the
        shared main memory, performs simple address translations.

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


class MemManager():
    def __init__(self, base_addr=0):
        self.__mmap__ = dict()
        self.__next_addr = base_addr

    def allocate(self, var, num_bytes):
        if var in self.__mmap__:
            print("mem  manager: warning duplicate allocation")

        self.__mmap__[var] = (self.__next_addr, num_bytes)
        self.__next_addr = self.__next_addr + num_bytes
        # print("next add: %d "%self.__next_addr)

    def deallocate(self, var):
        if not (var in self.__mmap__):
            print("mem manager: warning deallocation key does not exit")
            return

        del self.__mmap__[var]

    def get_address(self, var, k, size_t):

        if not (var in self.__mmap__):
            print("mem manager: warning requested key %s does not exist" % var)
            return -1

        add_begin = self.__mmap__[var][0]
        add_size = self.__mmap__[var][1]

        req_add = add_begin + k * size_t
        if (req_add > add_begin + add_size):
            print("mem manager: access out of bounds")

        return req_add
