# Using Mondrian Cross-Conformal Prediction to Predict Properties of Large Unbalanced Bioactivity Datasets
# Jiangming Sun, Lars Carlsson, Ernst Ahlberg, Ulf Norinder, Ola Engkvist and Hongming Chen, http://biorxiv.org/content/early/2017/03/16/116764

import os
import numpy as np
import numpy.ma as ma
from time import time
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from sklearn import cross_validation
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from numpy import linalg as LA
import scipy
import pandas as pd
np.random.seed()
linear = False

# Define all functions needed to do the various versions of CPs.

def zeroClassFirst(X,y):
	if (y[0]>0.5):
		for i,yLabel in enumerate(y[1:]):
			if (y[i]<0.5):
				yTmp = y[i]; XTmp = X[i]
				y[i] = y[0]; X[i] = X[0]
				y[0] = yTmp; X[0] = XTmp
				break

def computePValues(alpha_clf,X_properTrain,y_properTrain,X_calibration,y_calibration,X_test,outputdir):
	# Compute nonconformity scores
	# Do timings in here as well so that we can compare the different methods.
	trainTime = 0.0; predictTime = 0.0
	t0 = time()
	# The reason for doing this is that libsvm always uses the label of the first training 
	# example to define the negative side of the decision boundary, unless class labels are -1 and 1.
	zeroClassFirst(X_properTrain,y_properTrain)
	# Build model a calculate nonconformity scores for the calibration set. 
	alpha_clf.fit(X_properTrain,y_properTrain)
	y_calibrationAlphas = alpha_clf.decision_function(X_calibration)
	conditionZero = ma.masked_less_equal(y_calibration, 0.5)
	conditionOne = ma.masked_greater(y_calibration, 0.5)	
	if (y_properTrain[0] < 0.5 ): # The higher value response will have the higher decision value.
		alpha_zeros = np.extract(conditionZero.mask,y_calibrationAlphas)
		alpha_ones = np.extract(conditionOne.mask,-1.0*y_calibrationAlphas) # Negate to create a nonconformity score.
	else: # The lower value response will have the higher decision value.
		print("The first training example does not have a zero label!!!")
		sys.exit()
	alpha_zeros.sort()
	alpha_ones.sort()
	
	trainTime = trainTime + time() - t0
	# Compute p-values for the test examples.
	t0 = time()
	y_testAlphas = alpha_clf.decision_function(X_test)
	# Searching is done from the left, thus a larger value of searchsorted is more nonconforming.
	# Indexing start at 0 this is why we set +2 rather than +1.
	p_zeros = 1.0-1.0*(np.searchsorted(alpha_zeros,y_testAlphas)+1)/(len(alpha_zeros)+1)
	p_ones = 1.0-1.0*(np.searchsorted(alpha_ones,-1.0*y_testAlphas)+1)/(len(alpha_ones)+1)

	predictTime = predictTime + time() - t0
	
	# Print dots to indicate progress.
	print (".")

	return p_zeros,p_ones,trainTime,predictTime,y_testAlphas

def predictCCP(alpha_clf,X_train,y_train,X_test,outputdir):
	# The CCP is computing the mean of all p-values generated for the different folds for an example.
	totalTrainTime = 0.0; totalPredictTime = 0.0
	K = 5 # The number of folds. Se till att proper training och calibration motsvarar varandra mellan CCP, ACP!!!
	skf = cross_validation.StratifiedKFold(y_train, shuffle=True, n_folds=K, random_state=None)
	count = 0
	for properTrain_index, calibration_index in skf:
		X_properTrain, X_calibration = X_train[properTrain_index], X_train[calibration_index]
		y_properTrain, y_calibration = y_train[properTrain_index], y_train[calibration_index]
		p_0,p_1,trainTime,predictTime,alphatest = computePValues(alpha_clf,X_properTrain,y_properTrain,X_calibration,y_calibration,X_test,outputdir)
		if (count==0):
		  p_zeros = 1.0*p_0; p_ones = 1.0*p_1; at= 1.0*alphatest
		else:
		  p_zeros = np.vstack((p_zeros,p_0)); p_ones = np.vstack((p_ones,p_1)); at= np.vstack((at,alphatest))
		totalTrainTime = totalTrainTime + trainTime; totalPredictTime = totalPredictTime + predictTime
		count = count + 1
	return np.mean(p_zeros, axis=0), np.mean(p_ones, axis=0),totalTrainTime,totalPredictTime
	#return np.median(p_zeros, axis=0), np.median(p_ones, axis=0),totalTrainTime,totalPredictTime

def plotCVAccuracyGrid(grid_scores,C_vals,gamma_Vals,outputdir):
	# plot the scores of the grid
	# grid_scores_ contains parameter settings and scores
	# We extract just the scores
	len_gamma = len(gamma_Vals)
	len_c = len(C_vals)
	scores = [x[1] for x in grid_scores]
	scores = np.array(scores).reshape(len_c, len_gamma)
	gV = []
	for g in gamma_Vals:
		gV.append(expFormatter(g,None))
	cV = []
	for c in C_vals:
		cV.append(expFormatter(c,None))
	plt.figure(figsize=(8, 6))
	#fig, ax = plt.subplots()
	#ax.yaxis.set_major_formatter(FuncFormatter(roundFormatter))
	#ax.xaxis.set_major_formatter(FuncFormatter(expFormatter))
	plt.subplots_adjust(left=.2, right=0.95, bottom=0.15, top=0.95)
	plt.imshow(scores, interpolation='nearest', cmap=plt.cm.hot)
	plt.xlabel('gamma')
	plt.ylabel('C')
	plt.colorbar()
	plt.xticks(np.arange(len_gamma), gV, rotation=45)
	plt.yticks(np.arange(len(C_vals)), cV)
	plt.title('Validation accuracy')
	plt.savefig(outputdir+'/CVaccuracy.eps',format='eps')
	# plt.show()

def roundFormatter(x,pos):
	return '%.2f' % x

def expFormatter(x,pos):
	return '%.3g' % x

def plotCVAccuracyLine(grid_scores,C_vals):
	# plot the scores of the grid
	# grid_scores_ contains parameter settings and scores
	# We extract just the scores
	scores = [x[1] for x in grid_scores]
	scores = np.array(scores).reshape(len(C_vals))

	plt.figure(figsize=(8, 6))
	plt.subplots_adjust(left=.2, right=0.95, bottom=0.15, top=0.95)
	#plt.imshow(scores, interpolation='nearest', cmap=plt.cm.hot)
	plt.xlabel('C')
	plt.ylabel('Accuracy')
	plt.yticks(np.arange(len(scores)), scores)
	plt.xticks(np.arange(len(C_vals)), C_vals)
	plt.title('Validation accuracy')
	plt.savefig(outputdir+'/CVaccuracy.eps',format='eps')
#	plt.show()

def plotCalibration(p_0,p_1,y_t,filename):
	n = 100
	eps = 1.0*(np.arange(n)+1)/n
	error_rate = np.zeros(n)
	for k in range(0,n):
		m = len(y_t)
		for i in range(0,m):
			if (y_t[i] > 0.5): # One class true response
				if p_1[i] <= eps[k]: # One class not in prediction set
					error_rate[k] = error_rate[k] + 1.0/m                
			else: # Zero class true response
				if p_0[i] <= eps[k]: # Zero class not in prediction set
					error_rate[k] = error_rate[k] + 1.0/m
	plt.scatter(eps,error_rate)
	plt.axis([0,1,0,1])
	plt.xlabel("expected error")
	plt.ylabel("observed error")
	plt.title("Calibration plot")
	plt.plot([0, 1], [0, 1], 'k-')
	plt.savefig(filename,format='eps')
#	plt.show()
	return LA.norm(eps-error_rate,ord=None)

def returnConMatrix(p_0,p_1,y_t,alpha):
	tp=0
	fp=0
	tn=0
	fn=0
	uncertains=0
	emptys=0
	tp = sum(np.logical_and(np.logical_and(p_1>alpha, p_0<=alpha),y_t==1))
	fp = sum(np.logical_and(np.logical_and(p_1>alpha, p_0<=alpha),y_t==0))
	tn = sum(np.logical_and(np.logical_and(p_0>alpha, p_1<=alpha),y_t==0))
	fn = sum(np.logical_and(np.logical_and(p_0>alpha, p_1<=alpha),y_t==1))
	uncertains = sum(np.logical_and(p_0>alpha, p_1>alpha))
	emptys = sum(np.logical_and(p_0<=alpha, p_1<=alpha))
	return tp, fp, tn, fn, uncertains, emptys


import sys, getopt
def main(argv):
	# Overall timing
	t0 = time()
	# Let's read in a classification dataset in libsvm format.
	# from sklearn.datasets import load_svmlight_file
	inputfile = '' # traing set
	testfile = '' # test set
	outputdir = '' # output directory
	
	try:
		opts, args = getopt.getopt(argv,"hi:t:o:",["ifile=","testfile=","outputdir="])
	except getopt.GetoptError:
		print("mccp_openmp_sklearn.py -i <ifile> -t <testfile> -o <outputdir>")
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
		 print("mccp_openmp_sklearn.py -i <ifile> -t <testfile> -o <outputdir>")
		 sys.exit()
		elif opt in ("-i", "--ifile"):
		 inputfile = arg
		elif opt in ("-t", "--testfile"):
		 testfile = arg
		elif opt in ("-o", "--outputdir"):
		 outputdir = arg
	print(inputfile)
	print(testfile)
	print(outputdir)
	#create output dir
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	
	from sklearn.datasets import load_svmlight_file
	X, y = load_svmlight_file(inputfile)
	X_train = X
	y_train = y
	# Do a grid search based on cross validation to find optimum.

	if linear:
		C_vals = np.logspace(-2.0,1.0,10)
		tuned_parameters = [{'kernel': ['linear'], 'C': C_vals}]
	else:
		C_vals = np.logspace(0.0,10,num=6,base=2)
		gamma_vals = np.logspace(-1.0,-5.0,10)
		tuned_parameters = [{'kernel': ['rbf'], 'gamma': gamma_vals,'C': C_vals}]#,
	score = 'cohen_kappa_score'
	#score = 'precision'
	from sklearn.metrics import cohen_kappa_score, make_scorer
	kappa_scorer = make_scorer(cohen_kappa_score)

	print("# Tuning hyper-parameters for %s" % score)
	print()
	if linear:# Maste kanske explicit kalla pa liblinear varianten har.
		#clf = GridSearchCV(SVC(C=1,cache_size=4000,probability=True), tuned_parameters, cv=5, scoring='%s' % score, n_jobs=1)
		clf = GridSearchCV(SVC(C=1,cache_size=6400), tuned_parameters, cv=5, scoring=kappa_scorer, n_jobs=1)
	else:
		#clf = GridSearchCV(SVC(C=1,cache_size=4000,probability=True), tuned_parameters, cv=5, scoring='%s' % score, n_jobs=1)
		clf = GridSearchCV(SVC(C=1,cache_size=6400,class_weight='balanced'), tuned_parameters, cv=5, scoring=kappa_scorer, n_jobs=1)	
	clf.fit(X_train, y_train)
	print("Best parameters set found on training set:")
	print()
	print(clf.best_params_)
	print()
	print("Grid scores on training set:")
	print()
	for params, mean_score, scores in clf.grid_scores_:
		print("%0.3f (+/-%0.03f) for %r" % (mean_score, scores.std() * 2, params))
	print()
	print()

	#work on test data set
	X_test, y_test = load_svmlight_file(testfile)
	classifier = SVC(C=clf.best_params_['C'],kernel=clf.best_params_['kernel'],gamma=clf.best_params_['gamma'],class_weight='balanced',cache_size=6400)

	classifier_svm = SVC(C=clf.best_params_['C'],kernel=clf.best_params_['kernel'],gamma=clf.best_params_['gamma'],class_weight='balanced', cache_size=6400,probability=True)
	classifier_svm.fit(X_train,y_train)
	svmprob_values=classifier_svm.predict_proba(X_test)
	svm_predict_class=classifier_svm.predict(X_test) 

	print("CCP validity:"); p_0_CCP, p_1_CCP, trainTime_CCP, predictTime_CCP = predictCCP(classifier,X_train,y_train,X_test,outputdir)

	#ccp_rst=np.vstack((y_test,p_0_CCP,p_1_CCP)).reshape(3,len(y_test)).transpose()
	#active: 1 or inactive: 0
	active_label=np.argmax(np.vstack((p_0_CCP,p_1_CCP)), axis=0) 
	
	ccp_svm_rst=np.vstack((y_test,active_label,p_0_CCP,p_1_CCP,svm_predict_class,svmprob_values.transpose())).reshape(7,len(y_test)).transpose()
	fheader="Y_test_label\tMCCP_pred_label\tp_0_CCP\tp_1_CCP\tSvm_pred_label\tSvm_Prob_0\tSvm_Prob_1"
	np.savetxt(outputdir+'/mccp_svm_pred_rst.txt',ccp_svm_rst,fmt='%0.3f', delimiter='\t',header=fheader)

	plotCalibration(p_0_CCP, p_1_CCP, y_test,outputdir+'/ccpcalibration.eps')

	# Total execution time
	print("Total execution time:")
	print(time()-t0)

if __name__ == "__main__":
   main(sys.argv[1:])
