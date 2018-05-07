import numpy as np
import argparse
import random as rd
from math import sqrt
import matplotlib.pyplot as plt

# parser = argparse.ArgumentParser(description='Calculates the coefficients of ascend algorithm.', usage='python3 ascent_algorithm.py data_file')

# parser.add_argument('data_file', help='File of the experiments', type=str)
# parser.add_argument('--plot', help='Plot error graphs', type=str)
# args = parser.parse_args()




train_file = 'z3Vars.dat'
test_file = 'z3Vars.dat'
quasi = True
quasi_minimax = .05
stabilzation_factor = 1e-6
variables = None


plot = not True
map_train_data_X = None
map_train_data_Y = None
map_test_data_X = None
map_test_data_Y = None

variables = np.array([
[ 0, 0, 0 ],
[ 0, 0, 1 ],
[ 0, 1, 0 ],
[ 0, 1, 1 ],
[ 1, 0, 0 ],
[ 1, 0, 1 ],
[ 1, 1, 0 ],
[ 1, 1, 1 ]])
	


def get_approximant(B, f):
	#print(">> CALCULATING COEFFICIENTS...\n")
	C = np.dot(B,f)
	#print("--- CURRENT COEFFICIENTS FOUND: ")
	#print(C)
	return C,C[0] ## array containing the coefficients and the error [e, c1, c2, ..., cm]

## check wether the minimax norm is satisfied for the outter set (E)
def test_aproximant(outter_set, coefficients, outter_f): 
	#print(">> TESTING COEFFICIENTS IN THE OUTTER SET...\n")
	epsi_ph = -1000000000000000
	for i in range(outter_set.shape[0]):
		error = outter_f[i] - np.dot(outter_set[i,:], coefficients[1:]) ## excluding the error (coefficients[0])
		abs_err = abs(error)
		if abs_err > epsi_ph:
			epsi_ph = abs_err # error absoluto mayor
			IE = i # indice del error
			sgn = np.sign(error) # signo del error
	return epsi_ph,IE,sgn

## swap one element of the inner set with the outter set
def swap(outter_set, in_f, out_f, mu,IE, A,B):
	amp1 = np.insert(outter_set[IE], 0, mu, 0)
	lamb = np.dot(amp1,B)
	betha_max = -10000000000
	betha_i_max = -1
	for x in range(B.shape[1]):
		betha = mu*(lamb[x]/B[0,x])
		if betha > betha_max:
			betha_max = betha
			betha_i_max = x
	
	#print("\t SWAPPING: II = %d | EI = %d"%(betha_i_max, IE))
	tempo = np.array(A[betha_i_max,1:], copy=True)
	A[betha_i_max] = np.insert(outter_set[IE], 0, mu, 0)
	outter_set[IE] = tempo

	## swap in f too
	tempo = in_f[betha_i_max]
	in_f[betha_i_max] = out_f[IE]
	out_f[IE] = tempo

	B = get_new_inverse(B, betha_i_max, lamb)

	return A, B, outter_set, in_f, out_f

## by 4th method in the minimax theory document, the shape of inner set is m+1 X m
def get_signs_4(inner_set):
	#print(">> GETTING MINIMAX SIGNS...")
	inner_set = inner_set.T
	signs = np.linalg.solve(inner_set[:,:inner_set.shape[1]-1], inner_set[:,inner_set.shape[1]-1])
	#print("minimax signs by 4th method:")
	signs = np.insert(np.sign(signs), signs.shape[0], -1, 0)
	#print(signs)
	return signs

## map dataset into variables arrangement
def map(dataset, variables):
	#print(">> MAPPING DATA...\n")
	new_dataset = np.array([[]])
	for term in variables:
		new_col = np.ones(dataset.shape[0])
		i = 0 ## indice de la columna del dataset original
		for exp in term:
			new_col *= dataset[:,i]**exp
			i+=1
		new_dataset = np.append(new_dataset, [new_col])
	new_dataset = np.reshape(new_dataset,(variables.shape[0],dataset.shape[0])).T
	return new_dataset

## stabilize data
def stabilize(data, factor=0.000001):
	shape = data.shape
	return np.array([x+(rd.random()*stabilzation_factor) if (x==0) else x*(1+rd.random()*factor) for x in data.reshape(-1)]).reshape(shape)

## set order of data in relation to origin vector O
def set_order(data):
	origin_v = np.array([np.min(x) for x in data.T])
	distances = np.array([np.sum((x-origin_v)**2) for x in data])
	data = data[np.argsort(distances)]
	return data

## efficient way to calculate new inverse
def get_new_inverse(B, betha, lambs):
	for i in range(B.shape[0]):
		B[i,betha] = B[i,betha]/lambs[betha]
	for i in range(B.shape[0]):
		for j in range(B.shape[0]):
			if i!=betha:
				B[j,i] = B[j,i]-lambs[i]*B[j,betha]
	return B

## testing coefficients with test file
def test_rms(C):
	tst_rms = 0
	for x in range(test_data_X.shape[0]):
		tst_rms += (test_data_Y[x]-np.dot(C[1:].T,test_data_X[x]))**2
	return sqrt(tst_rms/test_data_X.shape[0])

## ascend algorithm
def ascend(m):

	data_ids = np.array(range(train_data_X.shape[0]))
	outter_ids = data_ids
	
	## eleccion aleatoria
	inner_ids = np.random.choice(data_ids, size=m+1, replace=False)
	## elegir los primeros m+1 elementos(para prueba)
	#inner_ids = [x for x in range(m+1)]
	## elegir verctores equidistantes despues de ordenar
	#inner_ids = np.linspace(0,train_data_X.shape[0]-1,m+1, dtype=int)	
	
	outter_ids = np.delete(outter_ids, inner_ids)

	#print(">> SPLITTING DATA INTO OUTTER AND INNER SET")
	inner_set = train_data_X[inner_ids,:]
	inner_f = train_data_Y[inner_ids,:]
	outter_set = train_data_X[outter_ids,:]
	outter_f = train_data_Y[outter_ids,:]

	signs = get_signs_4(inner_set)
	A = np.insert(inner_set, 0, signs, axis=1)

	#print(">> CALCULATING FIRST INVERSE OF A => B\n")
	B = np.linalg.inv(A)
	iteration = 1

	EpsTh = np.array([])
	EpsPh = np.array([])

	while True:
		coefficients,int_err = get_approximant(B, inner_f)
		ext_err,Ie,mu = test_aproximant(outter_set, coefficients, outter_f)
		
		print("IT[%d] EpsTh(in):%f \tEpsPh(out):%f  "%(iteration, int_err, ext_err))
		EpsTh = np.insert(EpsTh,EpsTh.shape[0],int_err)
		EpsPh = np.insert(EpsPh,EpsPh.shape[0],ext_err)
		if (int_err >= ext_err) or (quasi and abs(int_err-ext_err)<=quasi_minimax):
			#if(quasi):
				#print("Quasi Minimax condition satisfied")
			break
		else:
			A,B,outter_set,inner_f,outter_f = swap(outter_set, inner_f, outter_f,mu, Ie, A,B)
		iteration += 1
		#ok = input("\n>> Continue?(y/n) ")
		#if ok=="n": break
	print("IT[%d] EpsTh(in):%f \tEpsPh(out):%f  "%(iteration, int_err, ext_err))


	# train rms
	trn_rms = 0
	for x in range(train_data_X.shape[0]):
		trn_rms += (train_data_Y[x]-np.dot(coefficients[1:].T,train_data_X[x]))**2
	trn_rms = sqrt(trn_rms/train_data_X.shape[0])

	# test rms
	tst_rms = test_rms(coefficients)


	# print("@IT: %d"%(iteration))
	# print("TRAIN RMS error: %f"%(trn_rms))
	# print("TEST RMS error: %f"%(tst_rms))
	# print("Inner Error: %f\tOutter Error: %f"%(int_err,ext_err))
	
	if plot:
		epsth, = plt.plot(EpsTh, label="Epsilon Theta")
		epsph, = plt.plot(EpsPh, label="Epsilon Phi")
		plt.legend(handles=[epsth, epsph])
		plt.text(EpsTh.shape[0]/2, max(EpsPh)/2, r"$\epsilon_{\theta}=%f,\ \epsilon_{\phi}=%f$"%(int_err, ext_err))
		plt.xlabel('iteraciones')
		plt.ylabel('errores')
		plt.show()
	
	return coefficients, trn_rms,tst_rms

def set():
	global train_data_X
	global train_data_Y
	global test_data_X
	global test_data_Y

	data = np.loadtxt(train_file)
	#data = set_order(data)
	
	train_data_X = data[:,:data.shape[1]-1]
	train_data_X = data[:,:data.shape[1]-1]
	train_data_X = map(train_data_X, variables)
	train_data_Y = data[:,data.shape[1]-1:]
	
	#print(">> STABILIZING DATA...")
	train_data_X = stabilize(train_data_X)
	train_data_Y = stabilize(train_data_Y)

	data = np.loadtxt(test_file)
	test_data_X = data[:,:data.shape[1]-1]
	test_data_X = map(test_data_X, variables)
	test_data_Y = data[:,data.shape[1]-1:]

def run():
	coefficients, trn_rms,tst_rms = ascend(variables.shape[0])
	return coefficients,trn_rms,tst_rms



def go():
	coefficients, trn_rms,tst_rms = ascend(variables.shape[0])
	print("coefficients")
	print(coefficients)
	print("train rms:")
	print(trn_rms)
	print("test rms:")
	print(tst_rms)
	#data = np.loadtxt(args.data_file)
	## ejemplo tomado de 6_FLTS.TXT para M=5 terminos
	#variables = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
	#                      [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
	#                      [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0],
	#                      [0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
	#                      [4,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]])
	## ejemplo tomado de TRN.TXT, archivo de entrenamiento de DB10-iris, para M=3 terminos
	#variables = np.array([[0,0,1,0],
	#                     [0,0,0,1],
	#                     [0,1,0,0]])
	## ejemplo para z2variables.csv
	# variables = np.array([[0,0],
	#                      [0,1],
	#                      [0,2],
	#                      [1,0],
	#                      [1,1],
	#                      [1,2],
	#                      [2,0],
	#                      [2,1],
	#                      [2,2]])
	## ejemplo para z3Vars.dat como el del pdf

	## ejemplo para V4.txt con grados 1 1 1 1
	# variables = np.array([[0,0,0,0],
	#                       [0,0,0,1],
	#                       [0,0,1,0],
	#                       [0,0,1,1],
	#                       [0,1,0,0],
	#                       [0,1,0,1],
	#                       [0,1,1,0],
	#                       [0,1,1,1],
	#                       [1,0,0,0],
	#                       [1,0,0,1],
	#                       [1,0,1,0],
	#                       [1,0,1,1],
	#                       [1,1,0,0],
	#                       [1,1,0,1],
	#                       [1,1,1,0],
	#                       [1,1,1,1]])
	# coefficients, rms = ascend(data, variables.shape[0], variables)
	# print("Coeffs:")
	# print(coefficients)
	# return coefficients,rms

set()
go()
