from multiprocessing import Process
import subprocess
import time
from colorama import Fore, Back, Style
import os
import json


from sympy import *
from sympy.parsing import mathematica
# sa, x = symbols("a x")

# mathParse = mathematica.mathematica()
# print (type(mathParse))
# print (mathParse)
def genNumbers():
	'''
	'''
	pass

def workerFct(threadNb, debug):
	# subprocess.ca
	colorDict= {'0':Fore.RED, '1':Fore.GREEN, '2':Fore.BLUE, '3':Fore.YELLOW}
	if debug:
		print ("Debug On!!")

	FNULL = open(os.devnull, 'w')

	while True:
		print (colorDict[str(threadNb % 4) ] + "Mathematica Out from Thread number: ", threadNb)
		try:
			subprocess.call('./SO11_Masses.m ThreadNb-' + str(threadNb), shell = True, stdout = FNULL,  stderr=subprocess.STDOUT)
			time.sleep(0.1 + threadNb/10)

			#  Get fron out dict .
		except:
			break
			print ('No Bueno!!!!!!!!!!!!!!!!!!!!!!!!')

		print (Style.RESET_ALL)

if __name__=='__main__':

	processes = [Process(target=workerFct, args=(int(threadNb), False) ) for threadNb in range(16) ]


	dataDictOut = {'k':89000, 'zL':35}
	with open("dataIn.json", 'w') as outFile:
	    json.dump(dataDictOut, outFile)


	for proc in processes:
		proc.start()

	time.sleep(15)

	for proc in processes:
		proc.terminate()
