#!/usr/bin/env python

from data import *
from peak_picking import *
from plot_rovib_spectra import *
import json
from os import makedirs, listdir, path, getcwd as cwd
from matplotlib import *
from j_transitions import *
from models import *
from output_data import outputData as out
from custom_types import *
from pathlib import Path

# use if csv file in same dir as 'code' folder 
csv_path = Path(cwd()).parent.absolute()

# use if csv file in same dir as .py files
# csv_path = cwd()

# input name of rovib spectra file here (NO file extension!)
csv_name = 'dbr-toluene-spectra'
rovib_spectra_csv = rf'{csv_path}/{csv_name}.csv'

# default plot dpi
dpi = 600

# default plot size (width x height)
width = 14
height = 8
rc('figure', figsize=(width, height))


if __name__ == '__main__':
    
    # peaks_json_name = rf'example-analysis/data.json'
    # with open(peaks_json_name) as f:
    #     peaks = json.load(f) 
    # print(peaks)
    while True:
        project_name = input('\n\nPlease enter the name of the project: ')
        prompt = input(f'\nis this correct? [y/n]: ')
        if prompt == 'y':
            break
        elif prompt == 'n':
            continue
        else:
            print('\nexiting...')
            exit(1)
    
    # project directory
    i = 1
    dir_path = f'{csv_path}/{project_name}'
    for f in listdir(cwd()):
        if path.basename(f) == dir_path:
            if path.isdir(dir_path):
                dir_path = f'{dir_path}-({i})'
                i += 1
    makedirs(dir_path)

    df, x_name, y_name =  loadCsv(rovib_spectra_csv)
    
    print(
        f'\nx column name: {x_name}',
        f'\ny column name: {y_name}'
        )
    peaks = init()
    
    # set show_plot = True if you want to peak pick
    show = False
    peak_pick(peaks, df, show, x_name, y_name)
    
    # get transitions 
    j_transitions = calc_j_values(peaks)
    j_transitions.to_csv(
        f'{dir_path}/transitions.csv',
        index=False
    )
    
    
    # plot spectrum
    show = False # change to True to see plot
    rovib_spec = plot_rovib_spectra(df, peaks, show)
    basename  = 'rovib-spectra'
    extension = '.png'
    file_name = f'{dir_path}/{basename}{extension}'
    rovib_spec.savefig(file_name, dpi=600)
    
    
    # model energies
    for key in peaks.keys():
        
        
        # init analyzed data
        output = out(dir_path, peaks[key]['specie'])
        
        # append state variables
        state_vars = peaks[key]['state-variables']
        output.StateVars(4, state_vars)

        
        # zero order
        fig = calcModel(output, peaks, x_name, key, model_order='zero')
        basename = f'{key}-zero-order-fit'
        extension = '.png'
        file_name = f'{dir_path}/{basename}{extension}'
        fig.savefig(file_name, dpi=600)
        
        # first order
        fig = calcModel(output, peaks, x_name, key, model_order='first')
        basename = f'{key}-first-order-fit'
        extension = '.png'
        file_name = f'{dir_path}/{basename}{extension}'
        fig.savefig(file_name, dpi=600)

        # second order
        fig = calcModel(output, peaks, x_name, key, model_order='second')
        basename = f'{key}-second-order-fit'
        extension = '.png'
        file_name = f'{dir_path}/{basename}{extension}'
        fig.savefig(file_name, dpi=600)
        
        # output data 
        output.output()
        
    # print newline 
    print("\n")
