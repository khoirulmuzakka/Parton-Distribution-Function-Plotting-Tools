import ctypes as ct
import os
import numpy as np
import sys

pathLib = os.path.abspath("./libs/")
f123l = ct.CDLL(pathLib+'/libf123L.so') 



