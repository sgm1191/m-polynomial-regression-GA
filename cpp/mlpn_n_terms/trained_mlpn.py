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
input = np.array([1,0.161753964,0.089943331,0.04864474,0.040851372])
data = np.loadtxt("recall_data.csv",delimiter=',')
for record in data:
	X = np.insert(record,0,1)
	y_cap = scale(predict(X))
	print("%f %d"%(y_cap,np.round(y_cap)))

