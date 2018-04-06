import numpy as np
import ascent as asc
from math import ceil,sqrt
from random import random as rd
from random import randint as ri


asc.train_file = 'DB24-glass/TRAIN.TXT'
asc.test_file = 'DB24-glass/TEST.TXT'
asc.plot = False
asc.stabilzation_factor = 1e-6
#asc.train_file = 'z3Vars_trn.dat'
#asc.test_file = 'z3Vars_tst.dat'
pc = 1
pm = .05
generations = 100
N = 100
n_terms = 5

n_vars = 9
max_deg = 3

L = n_terms*n_vars*ceil(np.log2(max_deg))
b2m = round(L*N*pm)

def set():
	global L
	global b2m
	L = n_terms*n_vars*ceil(np.log2(max_deg))
	b2m = round(L*N*pm)

## genera una poblacion aleatoria de 1s y 0s con probabilidad uniforme
def gen_pop():
	pop = np.random.choice([0, 1], size=N*L, p=[.5, .5])
	return pop.reshape(N,L)

## decodifica un individuo a una matriz txv donde 
## t es el numero de terminos y v el numero de variables
def decode(individual):
	var_exp = np.array([])
	bits_exp = ceil(np.log2(max_deg)) # bits oara representar el grado maximo
	term_len = n_vars*bits_exp
	for t in range(n_terms):
		term = np.array([])
		c_ind = t*term_len
		term_str = np.array2string(individual[c_ind:c_ind+term_len], separator='')[1:-1]
		for v in range(n_vars):
			init = v*bits_exp
			dec_exp = int(term_str[init:init+bits_exp],2)
			term = np.insert(term,term.size,dec_exp)
		var_exp = np.insert(var_exp,var_exp.shape[0],term)
	return var_exp.reshape(n_terms,n_vars)

## evalua los coeficientes(sin el error incluido), con el conjunto de prueba
def test_rms(C, terms):
	return asc.test_rms(C, test_file, terms)

## evalua cada individuo de la poblacion
def evaluate(population):
	fitness = np.array([])
	for ind in population:
		terms = decode(ind)
		asc.variables = terms
		asc.set()
		C,trn_rms,tst_rms = asc.run()
		fitness = np.insert(fitness, fitness.size, tst_rms)
	return fitness

## Deterministic Selection Annular Crossover, para poblacion duplicada
def annular_cross(population):
	for i in range(N//2):
		p = rd()
		if p <= pc:
			p = ri(1, L//2)
			c = np.array(population[i,p:])
			population[i,p:] = population[N-i-1,p:]
			population[N-i-1,p:] = c
	return population


def mutate(population):
	for i in range(b2m):
		p1 = round(rd()*L)-1
		p2 = round(rd()*N)-1
		if population[p2,p1] == 0:
			population[p2,p1] = 1
		else:
			population[p2,p1] = 0
	return population


def run():
	asc.set()
	# create random population
	population = gen_pop()
	# evaluate population
	fitness = evaluate(population)
	for g in range(generations):
		# duplicate population (& fitness)
		population = np.tile(population[:N,:],(2,1))
		fitness = np.tile(fitness[:N],2)
		# cross
		population = annular_cross(population)
		# mutation
		population = mutate(population)
		# evaluate new population
		fitness[:N] = evaluate(population[:N,:])
		# sort population
		sorted_indices = np.argsort(fitness)
		fitness = fitness[sorted_indices]
		population = population[sorted_indices,:]
		# print results
		print("gen[%d] best fitness = %f\t worst fitness = %f"%(g+1,fitness[0],fitness[N]))
		#print(decode(population[0]))

	return decode(population[0]),fitness[0]

terms,best_fit = run()
asc.variables = terms
asc.set()
C,trn_rms,tst_rms = asc.run()
print("terms")
print(terms)
print("coeficientes")
print(C)
print("train rms = ")
print(trn_rms)
print("test rms = ")
print(tst_rms)
