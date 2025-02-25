from __future__ import print_function, division
import sys
import argparse


print ("WARNING: Please be aware that ERRORS can arise if relative paths are use in the main.py file (to point to the config file) or in the config.ini file (to point to other files used in the analysis).")

parser = argparse.ArgumentParser(
    description='''
    Main module for X-PSI Importance Sample.
    It requires a file main_IS.py in the same directory as the main.py where the importance sampling settings are set. For changes in the priors, this need to include importing a different CustomPrior, e.g. CustomPrior_IS.py in the main and the presence in both CustomPrior.py and the CustomPrior_IS.py of the method density, to evaluate prior density functions. 
    ''',
    fromfile_prefix_chars='@')

parser.add_argument('--name-original-main', type=str, help='Name of original main file')
parser.add_argument('--name-importance-sample-main', type=str, help='Name of importance sample main file')
parser.add_argument('--path-to-mains', type=str, help='Absolute or relative path to main files')
parser.add_argument('--path-to-outputs', type=str, help='Absolute or relative path to output files')
parser.add_argument('--output-root-name', type=str, help='Name of output root file, no extension')
parser.add_argument('--names', nargs='+',help='List of names of model variables')
parser.add_argument('--likelihood-change', type=str, help='Importance sample for chaning likelihood')
parser.add_argument('--prior-change', type=str, help='Importance sample for chaning priors')
parser.add_argument('--weight-threshold', type=float, help='Weight threshold for importance sampling')
parser.add_argument('--overwrite', type=str, help='Overwrite')

args = parser.parse_args()

if not args.path_to_mains:
    args.path_to_mains = input('Specify Absolute or relative path to main files: ')

if not args.name_original_main:
    args.name_original_main = input('Specify name of original main file: ')
    
if not args.name_importance_sample_main:
    args.name_importance_sample_main = input('Specify name of importance sample main file: ')

if not args.path_to_outputs:
    args.path_to_outputs = input('Specify Absolute or relative path to output files: ')
    if args.path_to_outputs[-1] =='/':
        args.path_to_outputs = args.path_to_outputs[:-1]
        
if not args.output_root_name:
    args.output_root_name = input('Specify name of output root file, without extension: ')

if not args.names:
    names = input('Specify the name of the model variables: ') #check if input
    args.names = names.split()
else:
    args.names = args.names[0].split()

def check_prior_change(test):
    if not args.prior_change:
        args.prior_change = input('Specify if importance sampling changes priors of model variables  (True/False): ')
        
    if args.prior_change.capitalize() == 'False':
        args.prior_change=False
        test = False
        return test
    elif args.prior_change.capitalize() == 'True':
        args.prior_change=True
        test = False
        return test
    else:
        print ('Only True or False options')
        args.prior_change = input('Specify if importance sampling changes priors of model variables  (True/False): ')
        
def check_likelihood_change(test):
    if not args.likelihood_change:
        args.likelihood_change = input('Specify if importance sampling changes likelihood of model variables  (True/False): ')
        
    if args.likelihood_change.capitalize() == 'False':
        args.likelihood_change=False
        test = False
        return test
    elif args.likelihood_change.capitalize() == 'True':
        args.likelihood_change=True
        test = False
        return test
    else:
        print ('Only True or False options')
        args.likelihood_change = input('Specify if importance sampling changes likelihood of model variables  (True/False): ')
        
def check_overwrite(test):
    if not args.overwrite:
        args.overwrite = input('Specify if to rewrite the output file(s)  (True/False): ')
        
    if args.overwrite.capitalize() == 'False':
        args.overwrite=False
        test = False
        return test
    elif args.overwrite.capitalize() == 'True':
        args.overwrite=True
        test = False
        return test
    else:
        print ('Only True or False options')
        args.overwrite = input('Specify if to rewrite the output file(s)  (True/False): ')

TEST = True
while TEST:
    TEST = check_prior_change(TEST)
    
TEST = True
while TEST:
    TEST = check_likelihood_change(TEST)
    
TEST = True
while TEST:
    TEST = check_overwrite(TEST)

if not args.weight_threshold:
    args.weight_threshold = 1.0e-30
    
    

import sys
sys.path.append(args.path_to_mains)
import importlib
original_main = importlib.import_module(args.name_original_main)
importance_sampling_main = importlib.import_module(args.name_importance_sample_main)

from xpsisilva.Sample import importance as importance_sample

importance_sample(importance_sampling_main.likelihood,
                  original_main.likelihood,
                  args.path_to_outputs + '/' + args.output_root_name,
                  names = args.names,
                  likelihood_change = args.likelihood_change,
                  prior_change = args.prior_change,
                  weight_threshold = args.weight_threshold,
                  overwrite = args.overwrite)
