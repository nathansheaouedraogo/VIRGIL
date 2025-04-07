import pandas as pd
from math import pi, sqrt
from scipy.constants import h, c, physical_constants as const

def loadCsv(
    file_path: str
    ) -> tuple[pd.DataFrame, str, str]:
        
    df = pd.read_csv(file_path)
    col_names = list(df)
    x_name = col_names[0]
    y_name = col_names[1]
    
    return df, x_name, y_name

def calcConstant(
    r_val: float,
    p_val : float
) -> tuple[float, float]:
    avg = (r_val + p_val)/2
    sd = sqrt(((p_val-avg)**2 + (r_val-avg)**2)/2)
    se = sd/sqrt(2)
    return avg, se

def appendCalculatedConstant(
    peaks: dict,
    key: str,
    constant: str,
    order: str,
    r_val: float, 
    p_val: float
):
    avg, se = calcConstant(r_val, p_val)
    peaks[key]['P-Branch'][constant][order].append(r_val)
    peaks[key]['R-Branch'][constant][order].append(p_val)
    peaks[key]['constants'][constant][order][0].append(avg)
    peaks[key]['constants'][constant][order][1].append(se)
    return avg, se

def reducedMass(m1: float, m2: float) -> float:
    mu = ((m1*m2) / (m1+m2)) * const['atomic mass constant'][0]
    return mu 
    
def init():
    
    peaks = {}
    i = 1
    while True:
        
        # get specie
        while True:
            specie = input(f'\nPlease enter specie {i}: ')
            print(f'\nSpecie: {specie}')
            prompt = input(f'\nis this correct? [y/n]: ')
            if prompt == 'y':
                print('\n')
                break
            elif prompt == 'n':
                continue
            else:
                print('\nexiting...')
                exit(1)
        
        # get k
        while True:
            k = float(input(f'\nPlease enter k (N/m) for {specie}: '))
            print(f'\nSpring constant (k): {k} N/m')
            prompt = input(f'\nis this correct? [y/n]: ')
            if prompt == 'y':
                print('\n')
                break
            elif prompt == 'n':
                continue
            else:
                print('\nexiting...')
                exit(1)

        # get mu
        while True:
            m1_amu = float(input(f'\nPlease enter the mass (amu) of molecule 1: '))
            m2_amu = float(input(f'\nPlease enter the mass (amu) of molecule 2: '))
            mu = reducedMass(m1_amu, m2_amu)
            print(f'\nReduced Mass (mu): {mu} kg')
            prompt = input(f'\nis this correct? [y/n]: ')
            if prompt == 'y':
                mu_amu = mu / const['atomic mass constant'][0]
                m1 = m1_amu * const['atomic mass constant'][0]
                m2 = m2_amu * const['atomic mass constant'][0]
                print('\n')
                break
            elif prompt == 'n':
                continue
            else:
                print('\nexiting...')
                exit(1)    

        # calc estimated vibrational fundamental 
        f_naught =  (1/(2*pi)) * sqrt(k/mu) * 10**(-12) # freq
        nu_naught = ((f_naught*10**(12)) / (c*100))  
        print(f'\nEstimated Fundamental Freq. of Vib:      {f_naught} THz')
        print(f'\nnEstimated Fundamental Wavenumber of Vib: {nu_naught} cm^(-1)')
        print('\n')

        # get colors
        while True:
            center_hex = input(f'Please input the graph colour for {specie} in hex: ')
            # Convert the center colour hex to an integer.
            center_int = int(center_hex, 16)
            # Calculate branch colour as an integer.

            # Extract RGB values from the center colour (convert each two-character hex to an int).
            r1 = int(center_hex[0:2], 16)
            g1 = int(center_hex[2:4], 16)
            b1 = int(center_hex[4:], 16)

            # Convert branch_int back to a 6-digit hex string.
            #branch_hex = format(branch_int, '06x')
            r2 = ((r1 & 0xf7) << 1) & 0xff
            g2 = ((g1 & 0xf7) << 1) & 0xff
            b2 = ((b1 & 0xf7) << 1) & 0xff
            branch_hex = '{:02x}{:02x}{:02x}'.format(r2, g2, b2)

            RESET = '\033[0m'
            # For example, using a white background (255,255,255); change as needed.
            fmt0 = '\033[48;2;255;255;255m'
            fmt1 = '\033[38;2;{};{};{}m'.format(r1, g1, b1)
            fmt2 = '\033[38;2;{};{};{}m'.format(r2, g2, b2)
            msg1 = 'This is the colour of the center'
            msg2 = 'This is the colour of the branches'
            out1 = '\n' + fmt0 + fmt1 + msg1 + RESET
            out2 = '\n' + fmt0 + fmt2 + msg2 + RESET 
            print(out1, out2)
            
            prompt = input('\nis this correct? [y/n]: ')
            if prompt == 'y':
                print('\n')
                break
            elif prompt == 'n':
                continue
            else:
                print('\nexiting...')
                exit(1)
        
        init = {
            
            'specie' : specie, 
            
            'state-variables' : {
                'spring-const' : [k, 'N·m⁻¹'], 
                'mass1-amu' : [m1_amu, 'amu'],
                'mass1-kg'  : [m1, 'kg'],
                'mass2-amu' : [m2_amu, 'amu'],
                'mass2-kg'  : [m2, 'kg'],
                'reduced-mass-amu' : [mu_amu, 'amu'],
                'reduced-mass-kg' : [mu, 'kg'],
                },
            
            'constants' : {
                'centrifugal-distortion-const' : {
                        'second-order-model' : ([],[]),
                    },
                'fNaughtTHz' : {
                        'estimated' : ([f_naught, []]),
                        'zero-order-model' : ([],[]),
                        'first-order-model' : ([],[]),
                        'second-order-model' : ([],[]),
                    },
                'moment-of-inertia-amu-angstroms' : {
                        'zero-order-model' : ([],[]),
                        'first-order-model' : ([],[]),
                        'second-order-model' : ([],[]),
                    },
                'nuTildeNaught' : {
                    'estimated' : ([nu_naught], []),
                    'zero-order-model' : ([],[]),
                    'first-order-model' : ([],[]),
                    'second-order-model' : ([],[]),
                },
                'radius-equilibrium-angstrom' : {
                        'zero-order-model' : ([],[]),
                        'first-order-model' : ([],[]),
                        'second-order-model' : ([],[]),
                    },
                'rovib-coupling-const' : {
                        'first-order-model' : ([],[]),
                        'second-order-model' : ([],[]),
                    },
                'rotational-constant' : {
                        'zero-order-model' : ([],[]),
                        'first-order-model' : ([],[]),
                        'second-order-model' : ([],[]),
                    },
                },
            
            'P-Branch' : {
                'v_tilde(cm^-1)' :  [],
                'abs' : [],
                'j_values': {
                    'j_init' : [],
                    'j_fini' : [],
                    },
                'centrifugal-distortion-const' : {
                        'second-order-model' : [],
                    },
                'fNaughtTHz' : {
                    'zero-order-model' : [],
                    'first-order-model' : [],
                    'second-order-model' : [],
                },
                'moment-of-inertia-amu-angstroms' : {
                        'zero-order-model' : [],
                        'first-order-model' : [],
                        'second-order-model' : [],
                },
                'nuTildeNaught' : {
                    'zero-order-model' : [],
                    'first-order-model' : [],
                    'second-order-model' : [],
                },
                'radius-equilibrium-angstrom' : {
                        'zero-order-model' : [],
                        'first-order-model' : [],
                        'second-order-model' : [],
                    },
                'rovib-coupling-const' : {
                        'first-order-model' :  [],
                        'second-order-model' :  [],
                    },
                'rotational-constant' : {
                        'zero-order-model' :  [],
                        'first-order-model' :  [],
                        'second-order-model' :  [],
                    },
                },
            
            'R-Branch' : {
                'v_tilde(cm^-1)' :  [],
                'abs' : [],
                'j_values': {
                    'j_init' : [],
                    'j_fini' : [],
                    },
                'centrifugal-distortion-const' : {
                        'second-order-model' : [],
                    },
                'fNaughtTHz' : {
                    'zero-order-model' : [],
                    'first-order-model' : [],
                    'second-order-model' : [],
                },
                'moment-of-inertia-amu-angstroms' : {
                        'zero-order-model' : [],
                        'first-order-model' : [],
                        'second-order-model' : [],
                },
                'nuTildeNaught' : {
                    'zero-order-model' : [],
                    'first-order-model' : [],
                    'second-order-model' : [],
                },
                'radius-equilibrium-angstrom' : {
                        'zero-order-model' : [],
                        'first-order-model' : [],
                        'second-order-model' : [],
                    },
                'rovib-coupling-const' : {
                        'first-order-model' :  [],
                        'second-order-model' :  [],
                    },
                'rotational-constant' : {
                        'zero-order-model' :  [],
                        'first-order-model' :  [],
                        'second-order-model' :  [],
                    },
                },
            
            'colour' : {
                'center' : f'#{center_hex}',
                'branch' : f'#{branch_hex}'
            }
        }
        
        # append to peaks dict
        key_name = specie.lower()
        peaks[key_name] = init
        
        # prompt user to keep going or exit
        prompt = input(f'\nAdd additional species? [y/N]: ')
        if prompt == 'y':
            print('\n')
            continue
        else:
            break
        
    return peaks
