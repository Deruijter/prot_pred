import pandas as pd
import numpy as np
import pylab as p
import csv as csv
import datetime
import warnings
import glob
import os
import sys
import collections
from shutil import copyfile
from optparse import OptionParser, OptionGroup
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn import svm, linear_model
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import *
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import RFECV


file_path_ids = '../data/training_ids.csv'
file_path_data_sets = '../feature_sets/'
out = '3_machine_learn_data_sets.out'

ml_fs_results = []

regressor_types = {}
rfr_regressors = {}
rir_regressors = {}
svm_regressors = {}
knn_regressors = {}
nnr_regressors = {}

rfr_regressors['rfr_default'] = RandomForestRegressor()	
rfr_regressors['rfr_ne25'] = RandomForestRegressor(n_estimators=25)
rfr_regressors['rfr_ne50'] = RandomForestRegressor(n_estimators=50)
rfr_regressors['rfr_ne100'] = RandomForestRegressor(n_estimators=100)
rfr_regressors['rfr_mf1'] = RandomForestRegressor(max_features=1)
rfr_regressors['rfr_mf2'] = RandomForestRegressor(max_features=2)
rfr_regressors['rfr_mf3'] = RandomForestRegressor(max_features=3)
rfr_regressors['rfr_mfSqrt'] = RandomForestRegressor(max_features='sqrt')
rfr_regressors['rfr_mfLog2'] = RandomForestRegressor(max_features='log2')
rfr_regressors['rfr_mss3'] = RandomForestRegressor(min_samples_split=3)
rfr_regressors['rfr_mss5'] = RandomForestRegressor(min_samples_split=5)
rfr_regressors['rfr_mss10'] = RandomForestRegressor(min_samples_split=10)
rfr_regressors['rfr_mss0.1'] = RandomForestRegressor(min_samples_split=0.1)
rfr_regressors['rfr_mss0.5'] = RandomForestRegressor(min_samples_split=0.5)
rfr_regressors['rfr_mss1.0'] = RandomForestRegressor(min_samples_split=1.0)
rfr_regressors['rfr_msl3'] = RandomForestRegressor(min_samples_leaf=3)
rfr_regressors['rfr_msl5'] = RandomForestRegressor(min_samples_leaf=5)
rfr_regressors['rfr_ne100_msl5'] = RandomForestRegressor(n_estimators=100, min_samples_leaf=5)
rfr_regressors['rfr_msl10'] = RandomForestRegressor(min_samples_leaf=10)
rfr_regressors['rfr_msl0.01'] = RandomForestRegressor(min_samples_leaf=0.01)
rfr_regressors['rfr_msl0.05'] = RandomForestRegressor(min_samples_leaf=0.05)
rfr_regressors['rfr_msl0.10'] = RandomForestRegressor(min_samples_leaf=0.1)
regressor_types['rfr'] = rfr_regressors

svm_regressors['svm_default'] = svm.LinearSVR()
svm_regressors['svm_epsilon0.0'] = svm.LinearSVR(epsilon=0)
svm_regressors['svm_epsilon0.2'] = svm.LinearSVR(epsilon=0.2)
svm_regressors['svm_epsilon0.3'] = svm.LinearSVR(epsilon=0.3)
svm_regressors['svm_epsilon0.4'] = svm.LinearSVR(epsilon=0.4)
svm_regressors['svm_epsilon0.5'] = svm.LinearSVR(epsilon=0.5)
svm_regressors['svm_epsilon0.6'] = svm.LinearSVR(epsilon=0.6)
svm_regressors['svm_epsilon0.7'] = svm.LinearSVR(epsilon=0.7)
svm_regressors['svm_epsilon0.8'] = svm.LinearSVR(epsilon=0.8)
svm_regressors['svm_epsilon0.9'] = svm.LinearSVR(epsilon=0.9)
svm_regressors['svm_epsilon1.0'] = svm.LinearSVR(epsilon=1.0)
svm_regressors['svm_c0.5'] = svm.LinearSVR(C=0.5)
svm_regressors['svm_c2.0'] = svm.LinearSVR(C=2.0)
svm_regressors['svm_lossSEI'] = svm.LinearSVR(loss='squared_epsilon_insensitive')
regressor_types['svm'] = svm_regressors

rir_regressors['rir_default'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='auto')
rir_regressors['rir_alph2.0'] = linear_model.Ridge(alpha=2.0, normalize=False,copy_X=True,solver='auto')
rir_regressors['rir_alph0.5'] = linear_model.Ridge(alpha=0.5, normalize=False,copy_X=True,solver='auto')
rir_regressors['rir_copyXfalse'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=False,solver='auto')
rir_regressors['rir_normTrue'] = linear_model.Ridge(alpha=1.0, normalize=True,copy_X=True,solver='auto')
rir_regressors['rir_slvsvd'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='svd')
rir_regressors['rir_slvcho'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='cholesky')
rir_regressors['rir_slvlsqr'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='lsqr')
rir_regressors['rir_slvscg'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='sparse_cg')
rir_regressors['rir_slvsag'] = linear_model.Ridge(alpha=1.0, normalize=False,copy_X=True,solver='sag')
regressor_types['rir'] = rir_regressors

knn_regressors['knn_default'] = KNeighborsRegressor(n_neighbors=5, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn003'] = KNeighborsRegressor(n_neighbors=3, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn010'] = KNeighborsRegressor(n_neighbors=10, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn010_wtd'] = KNeighborsRegressor(n_neighbors=10, weights='distance',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn050'] = KNeighborsRegressor(n_neighbors=50, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn150'] = KNeighborsRegressor(n_neighbors=150, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn500'] = KNeighborsRegressor(n_neighbors=500, weights='uniform',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn150_wtd'] = KNeighborsRegressor(n_neighbors=150, weights='distance',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
knn_regressors['knn_nn500_wtd'] = KNeighborsRegressor(n_neighbors=500, weights='distance',algorithm='auto',leaf_size=30,p=2,metric='minkowski',metric_params=None)
regressor_types['knn'] = knn_regressors

nnr_regressors['nnr_default'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='adam',max_iter=200)
nnr_regressors['nnr_hls150'] = MLPRegressor(hidden_layer_sizes=150,activation='relu',solver='adam',max_iter=200)
nnr_regressors['nnr_hls150_mi400'] = MLPRegressor(hidden_layer_sizes=150,activation='relu',solver='adam',max_iter=400)
nnr_regressors['nnr_hls300'] = MLPRegressor(hidden_layer_sizes=300,activation='relu',solver='adam',max_iter=200)
nnr_regressors['nnr_hls300_mi400'] = MLPRegressor(hidden_layer_sizes=300,activation='relu',solver='adam',max_iter=400)
nnr_regressors['nnr_hls500'] = MLPRegressor(hidden_layer_sizes=500,activation='relu',solver='adam',max_iter=200)
nnr_regressors['nnr_hls500_mi400'] = MLPRegressor(hidden_layer_sizes=500,activation='relu',solver='adam',max_iter=400)
nnr_regressors['nnr_hls1000'] = MLPRegressor(hidden_layer_sizes=1000,activation='relu',solver='adam',max_iter=200)
nnr_regressors['nnr_hls1000_mi400'] = MLPRegressor(hidden_layer_sizes=1000,activation='relu',solver='adam',max_iter=400)
nnr_regressors['nnr_mi0300'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='adam',max_iter=300)
nnr_regressors['nnr_mi0500'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='adam',max_iter=500)
nnr_regressors['nnr_mi1000'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='adam',max_iter=1000)
nnr_regressors['nnr_acttanh'] = MLPRegressor(hidden_layer_sizes=100,activation='tanh',solver='adam',max_iter=200)
nnr_regressors['nnr_actlogi'] = MLPRegressor(hidden_layer_sizes=100,activation='logistic',solver='adam',max_iter=200)
nnr_regressors['nnr_actiden'] = MLPRegressor(hidden_layer_sizes=100,activation='identity',solver='adam',max_iter=200)
nnr_regressors['nnr_slvlbfgs'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='lbfgs',max_iter=200)
nnr_regressors['nnr_slvsgd'] = MLPRegressor(hidden_layer_sizes=100,activation='relu',solver='sgd',max_iter=200)
regressor_types['nnr'] = nnr_regressors




feature_selectors = dict()
feature_selectors['sfm_en0.25_0.05'] = linear_model.ElasticNet(l1_ratio=0.25, alpha=0.05)
feature_selectors['sfm_en0.70_0.05'] = linear_model.ElasticNet(l1_ratio=0.70, alpha=0.05)
feature_selectors['sfm_en0.95_0.05'] = linear_model.ElasticNet(l1_ratio=0.95, alpha=0.05)
feature_selectors['sfm_et'] = ExtraTreesRegressor()
feature_selectors['sfm_ls'] = svm.LinearSVR()
feature_selectors['var_vt0.05'] = VarianceThreshold(0.05)
feature_selectors['rfe_et'] = RFECV(ExtraTreesRegressor(), step=0.05, n_jobs=2, cv=4, scoring='r2')
feature_selectors['rfe_ls'] = RFECV(svm.LinearSVR(), step=0.05, n_jobs=2, cv=4, scoring='r2')
feature_selectors['rfe_rr'] = RFECV(linear_model.Ridge(), step=0.05, n_jobs=2, cv=4, scoring='r2')

def output(out_text):
	dt = '{0:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
	out_text = '{0} - {1}'.format(dt, out_text)
	print out_text
	with open(out, 'a+') as outfile:
		outfile.write('{0}\n'.format(out_text))

def parse_commandline():
    usage = ''
    version = ''
    description = ''
    epilog = ''
    parser = OptionParser(usage=usage, description=description,
                          version=version, epilog=epilog)

    # sequence/alignment options:    
    parser.add_option("-d", "",  dest="data_set_file_name", type="string",
                     help="File name of data set")
    parser.set_defaults(data_set_file_name='')
    
    # get the options:
    (options, args) = parser.parse_args()

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))

    if not options.data_set_file_name:
    	output('Please pass data set file name using -d [filename]')
    	exit()
    
    # clean up (recommended):
    del(parser)
    return options


def run_clf2(clf, X, Y):
	score = cross_val_score(clf, X, Y, scoring='r2', cv=8, n_jobs=1).mean()
	return score

def run_clf2_final(clf, X, Y, X_test, Y_test, regressor_type):
	# Create model
	clf.fit(X,Y)

	feature_importance = ''
	if regressor_type == 'rir' or regressor_type == 'svm':
		feature_importance = clf.coef_ * -1
	elif regressor_type == 'rfr':
		feature_importance = clf.feature_importances_

	# Apply and score on test data
	return([clf.score(X_test,Y_test), feature_importance , clf.predict(X_test) ])

def create_feature_selection_sets():
	dataset_paths = glob.glob("{0}data_set_ultra_mega2.csv".format(file_path_data_sets))

	# read output file (in case we want to check things or not run the same thing twice)
	ml_results = dict()
	with open('../ml_output/ml_scores.csv', 'r+') as csv_file:
		reader = csv.reader(csv_file)
		for row in reader:
			ml_results['{0},{1},{2}'.format(row[0],row[1],row[2])] = '{0},{1},{2}'.format(row[3],row[4],row[5])

	fs_features = dict()
	with open('../ml_output/fs_features.csv', 'r+') as csv_file:
		reader = csv.reader(csv_file)
		for row in reader:
			fs_features['{0},{1}'.format(row[0],row[1])] = '{0}'.format(row[2])

	for dataset_path in dataset_paths:
		data = pd.read_csv(dataset_path, sep='\t', header=0)
		ids = pd.read_csv(file_path_ids, sep='\t', header=0)
		data_set_name = os.path.basename(dataset_path)
		data_set_name = data_set_name[9:-4]

		for i in range(1,1+10):
			print "{0:%Y-%m-%d %H:%M:%S} -  ITERATION {1}".format(datetime.datetime.now(),i)
			ids_train = ids[ids['set_{0}'.format(i)] == True]
			data_train = data[data['sys_id'].isin(ids_train['sys_id'])]

			ids_test = ids[ids['set_{0}'.format(i)] == False]
			data_test = data[data['sys_id'].isin(ids_test['sys_id'])]

			data_train = data_train.drop('sys_id', axis=1)
			data_test = data_test.drop('sys_id', axis=1)
			scaler = StandardScaler()
			scaler.fit(data_train)
			data_train = scaler.transform(data_train)
			data_train = pd.DataFrame(data_train)
			data_test = scaler.transform(data_test)
			data_test = pd.DataFrame(data_test)

			for fs_name, fs in feature_selectors.iteritems():

				fs_type = fs_name[0:3]
				data_set_name2 = '{0}_{1}'.format(data_set_name, fs_name)
				print "{0:%Y-%m-%d %H:%M:%S} -  {1}".format(datetime.datetime.now(),fs_name)

				X = data_train.ix[:,1:]
				Y = data_train.ix[:,0]
				X_test = data_test.ix[:,1:]
				Y_test = data_test.ix[:,0]

				features = ''

				if fs_type == 'sfm': # SELECT FROM MODEL [elastic net, extra trees, linear_svr]
					clf = fs.fit(X,Y)
					model = SelectFromModel(clf, prefit=True)
					X = model.transform(X)
					X_test = model.transform(X_test)
					features = model.get_support()

				elif fs_type == 'var': # VARIANCE THRESHOLD
					model = fs.fit(X,Y)
					X = fs.transform(X)
					X_test = fs.transform(X_test)
					features = fs.get_support()

				elif fs_type == 'rfe': # RFE [extra trees, linear_svr, ridge]
					model = fs.fit(X,Y)
					X = fs.transform(X)
					X_test = fs.transform(X_test)
					features = fs.get_support()

				fs_features['{0},{1}'.format(data_set_name2,i)] = features
				fs_features_ordered = collections.OrderedDict(sorted(fs_features.items(), key=lambda t: t[0]))
				with open('../ml_output/fs_features.csv', 'r+') as csv_file:
					for key, value in fs_features_ordered.items():
						csv_file.write('{0},{1}\n'.format(key, value))

				ml_results = build_learn_validate(X,Y,X_test,Y_test,data_set_name, ml_results)

			# output to file
			ml_results_ordered = collections.OrderedDict(sorted(ml_results.items(), key=lambda t: t[0]))
			with open('../ml_output/ml_scores.csv', 'r+') as csv_file:
				for key, value in ml_results_ordered.items():
					csv_file.write('{0},{1}\n'.format(key, value))

def build_learn_validate(X,Y,X_test,Y_test,data_set_name, i, ml_results):
	for regressor_type, ml_regressors in regressor_types.items():
		best_learner = ''
		training_score_best = -999
		training_score = -999	
		regressors_done = 0	

		for regressor_name, ml_regressor in ml_regressors.iteritems():

			key = '{0},{1},{2}'.format(regressor_name,data_set_name,i)
			if key in ml_results:
				#output("SKIP {0} - {1} -> already done".format(key, i))
				vals = ml_results[key].split(',')
				training_score = vals[0]
				test_score = vals[1]
				feature_importance = vals[2]
				regressors_done = regressors_done + 1

			else:
				#output("{0} on {2} - {1}".format(regressor_name, i, data_set_name))
				try:
					training_score = run_clf2(ml_regressor, X, Y)
					test_score = -999
					feature_importance = ''
				except:
					output("ERROR: {0}".format(sys.exc_info()[0]))
					training_score = -999
					test_score = -999
					feature_importance = ''

			if training_score > training_score_best:
				training_score_best = training_score
				best_learner = regressor_name

			ml_results[key] = '{0},{1},{2}'.format(training_score, test_score, feature_importance)

		# If not yet previously calculated
		if regressors_done < len(ml_regressors):
			#output("GET TEST SCORES {0}".format(best_learner))
			test_results = run_clf2_final(ml_regressors[best_learner], X, Y, X_test, Y_test, regressor_type)
			test_score = test_results[0]
			feature_importance = test_results[1]

			# redefine best record to include test score and coefs
			key = '{0},{1},{2}'.format(best_learner,data_set_name,i)
			ml_results[key] = '{0},{1},{2}'.format(training_score_best, test_score, feature_importance)

			with open('../ml_output/predictions/ml_predictions_{0}_{1}_{2}.csv'.format(regressor_type, data_set_name, i), 'w+') as csv_file:
				for j in range(0,len(Y_test)-1):
					csv_file.write('{0},{1}\n'.format(Y_test[j], test_results[2][j]))

	return(ml_results)

def main():

	#warnings.filterwarnings("ignore", category=DeprecationWarning)
	#output("STARTED SCRIPT")

	#create_feature_selection_sets()

	#output("STARTING REGULAR SETS")

	options = parse_commandline()
	data_set_file_name = options.data_set_file_name

	#output("Running learners for data set {0}".format(data_set_file_name))

	dataset_path = "{0}{1}".format(file_path_data_sets, data_set_file_name)

	# read output file (in case we want to check things or not run the same thing twice)
	ml_results = dict()

	data = pd.read_csv(dataset_path, sep='\t', header=0)
	ids = pd.read_csv(file_path_ids, sep='\t', header=0)
	data_set_name = os.path.basename(dataset_path)
	data_set_name = data_set_name[9:-4]

	for i in range(1,1+10):
		output("{0} IT {1}".format(data_set_name, i))
		ids_train = ids[ids['set_{0}'.format(i)] == True]
		data_train = data[data['sys_id'].isin(ids_train['sys_id'])]

		ids_test = ids[ids['set_{0}'.format(i)] == False]
		data_test = data[data['sys_id'].isin(ids_test['sys_id'])]

		data_train = data_train.drop('sys_id', axis=1)
		data_test = data_test.drop('sys_id', axis=1)

		# SCALE
		scaler = StandardScaler()
		scaler.fit(data_train)
		data_train = scaler.transform(data_train)
		data_train = pd.DataFrame(data_train)
		data_test = scaler.transform(data_test)
		data_test = pd.DataFrame(data_test)

		X = data_train.ix[:,1:]
		Y = data_train.ix[:,0]
		X_test = data_test.ix[:,1:]
		Y_test = data_test.ix[:,0]

		ml_results = build_learn_validate(X, Y, X_test, Y_test, data_set_name, i, ml_results)

		# output to file
		ml_results_ordered = collections.OrderedDict(sorted(ml_results.items(), key=lambda t: t[0]))

	with open('../ml_output/scores/ml_scores_{0}.csv'.format(data_set_name), 'w+') as csv_file:
		for key, value in ml_results_ordered.items():
			csv_file.write('{0},{1}\n'.format(key, value))
	output("{0} DONE! :)".format(data_set_name))


if __name__ == "__main__":
    main()
