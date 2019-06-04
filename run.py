from multiprocessing import Pool
import os
import numpy as np
import time


def drive(i):
	 os.system(f'./nematic {i}')

if __name__ == '__main__':
	start_time = time.time()
	temps = np.linspace(0.01,1,12)
	pool = Pool()
	pool.map(drive, temps)
	pool.close()
	pool.join()
	print("--- %s seconds ---" % (time.time() - start_time))