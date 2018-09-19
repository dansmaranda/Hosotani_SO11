from multiprocessing import Process
import subprocess
import time
from colorama import Fore, Back, Style
import os


def workerFct(threadNb, debug):
	# subprocess.ca
	colorDict= {'0':Fore.RED, '1':Fore.GREEN, '2':Fore.BLUE, '3':Fore.YELLOW}
	if debug:
		print ("Debug On!!")
	
	FNULL = open(os.devnull, 'w')

	while True:
		print (colorDict[str(threadNb % 4) ] + "Mathematica Out from Thread number: ", threadNb)
		try:
			subprocess.call('wolframscript -script SO11_Masses.m', shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
			time.sleep(0.1 + threadNb/10)
		except:
			break
			print ('No Bueno!!!!!!!!!!!!!!!!!!!!!!!!')

		print (Style.RESET_ALL)

if __name__=='__main__':
	processes = [Process(target=workerFct, args=(int(threadNb), False) ) for threadNb in range(16) ]

	for proc in processes:
		proc.start()

	time.sleep(15)

	for proc in processes:
		proc.terminate()
