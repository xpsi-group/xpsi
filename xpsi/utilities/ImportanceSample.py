from __future__ import print_function, division
import sys
import argparse

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
parser.add_argument('--likelihood-change', type=bool, help='Importance sample for chaning likelihood')
parser.add_argument('--prior-change', type=bool, help='Importance sample for chaning priors')
parser.add_argument('--weight-threshold', type=float, help='Weight threshold for importance sampling')
parser.add_argument('--overwrite', type=bool, help='Overwrite')

args = parser.parse_args()

if not args.path_to_mains:
    args.path_to_mains = raw_input('Specify Absolute or relative path to main files: ')

if not args.name_original_main:
    args.name_original_main = raw_input('Specify name of original main file: ')
    
if not args.name_importance_sample_main:
    args.name_importance_sample_main = raw_input('Specify name of importance sample main file: ')

if not args.path_to_outputs:
    args.path_to_outputs = raw_input('Specify Absolute or relative path to output files: ')
    if args.path_to_outputs[-1] =='/':
        args.path_to_outputs = args.path_to_outputs[:-1]
        
if not args.output_root_name:
    args.output_root_name = raw_input('Specify name of output root file, without extension: ')

if not args.names:
    names = input('Specify the name of the model variables: ') #check if input
    args.names = names.split()
else:
    args.names = args.names[0].split()

if not args.likelihood_change:
    while True:
        input_likelihood_change = raw_input('Specify if importanca sampling changes likelihood (True/False): ')
        if input_likelihood_change.capitalize() == 'False':
            args.likelihood_change=False
            break
        elif input_likelihood_change.capitalize() == 'True':
            args.likelihood_change=True
            break
        else:
            print ('Only True or False options')
        
if not args.prior_change:
    while True:
        input_prior_change = raw_input('Specify if importanca sampling changes priors of model variables  (True/False): ')
        if input_prior_change.capitalize() == 'False':
            args.input_prior_change=False
            break
        elif input_prior_change.capitalize() == 'True':
            args.input_prior_change=True
            break
        else:
            print ('Only True or False options')

if not args.weight_threshold:
    args.weight_threshold = 1.0e-30
    
if not args.overwrite:
    args.overwrite = True
    

import sys
sys.path.append(args.path_to_mains)
import importlib
original_main = importlib.import_module(args.name_original_main)
importance_sampling_main = importlib.import_module(args.name_importance_sample_main)

from xpsi.Sample import importance as importance_sample

importance_sample(importance_sampling_main.likelihood,
                  original_main.likelihood,
                  args.path_to_outputs + '/' + args.output_root_name,
                  names = args.names,
                  likelihood_change = args.likelihood_change,
                  prior_change = args.prior_change,
                  weight_threshold = args.weight_threshold,
                  overwrite = args.overwrite)
