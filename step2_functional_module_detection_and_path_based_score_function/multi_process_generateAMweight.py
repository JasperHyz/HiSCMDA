import sys
import getopt
import functional_module_detection
from concurrent.futures import ProcessPoolExecutor


if __name__=='__main__':
    PRIMES = [x for x in range(1, 6)]
#     try:
#         opts, args = getopt.getopt(sys.argv[1:], "n:", ['number='])
#     except getopt.GetoptError:
#         print('global_LOOCV.py -n <number>')
#         sys.exit(2)
#     for opt, arg in opts:
#         if opt in ("-n", "--number"):
#             N = int(arg)

    with ProcessPoolExecutor() as pool:
        pool.map(functional_module_detection.func_, PRIMES)

