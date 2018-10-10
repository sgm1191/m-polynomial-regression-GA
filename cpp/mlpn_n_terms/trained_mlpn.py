import numpy as np

def sigmoid(x):
	return (1 / (1 + np.exp(-x)))

def predict(input):
	l1 = np.array([
		[9.93878,-1.277238],
		[12.794103,2.149855],
		[-32.754846,9.336392],
		[1.290772,-8.214893],
		[8.711606,-21.087847]
		])
	lo = np.array([-.018833,.452167,.531239])

	h1 = sigmoid(np.dot(input, l1))
	h1 = np.insert(h1,0,1)

	h2 = np.dot(h1,lo)

	return h2

def scale(x):
	return (x*(13.0207319789-1.0063476744))+1.0063476744

#input = np.array([1,0.247738208,0.09904867,0.023874229,0.99999939])
#input = np.array([1,4,118,8926,2.55])

# the input must have the following forder
# num. of attributes, num. of rows, information, compression ratio
## for example
# 0.247738208,0.09904867,0.023874229,0.99999939

data = np.loadtxt("no_sc_new_experiments_recall.csv",delimiter=',')
i = 0
for record in data:
	i += 1
	#to scale data into [0,1)
	#atributes
	record[0] = (record[0]-1.356385937)/(36.25644988969-1.356385937)
	#number of records
	record[1] = (record[1]-59.4135663696)/(3134.5324164552-59.4135663696)
	#efective information calculated by gpm4 compression algorithm
	record[2] = (record[2]+5305.1795614531)/(258920.293929573+5305.1795614531)
	#compresion ratio
	record[3] = (record[3]-2.2681166805)/(38.0070218056-2.2681166805)
	#add bias value at the begining
	X = np.insert(record,0,1)
	y_cap = scale(predict(X))
	print("calculated number of terms= %f, rounded to %d"%(y_cap,np.round(y_cap)))

# for one input
# atributes
#input[1] = (input[1]-1.356385937)/(36.25644988969-1.356385937)
# number of inputs
#input[2] = (input[2]-59.4135663696)/(3134.5324164552-59.4135663696)
# efective information calculated by gpm4 compression algorithm
#input[3] = (input[3]+5305.1795614531)/(258920.293929573+5305.1795614531)
# compresion ratio
#input[4] = (input[4]-2.2681166805)/(38.0070218056-2.2681166805)
#y_cap = scale(predict(input))
#print("calculated number of terms= %f, rounded to %d"%(y_cap,np.round(y_cap)))

