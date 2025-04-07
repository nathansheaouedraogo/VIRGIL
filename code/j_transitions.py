from data import *

def calc_j_values(peaks: dict) -> pd.DataFrame:
    
    output = {
        'specie' : [],
        'branch' : [],
        'wavenumber' : [],
        'abs' : [],
        'j_init' : [],
        'j_fini' : [],
    }
    
    for key in peaks.keys():
        
        # fill R-Branch transitions
        j = 0
        for i in range(len(peaks[key]['R-Branch']['abs'])):
            peaks[key]['R-Branch']['j_values']['j_init'].append(j)
            peaks[key]['R-Branch']['j_values']['j_fini'].append(j+1)
            output['branch'].append('R-Branch')
            output['j_init'].append(j)
            output['j_fini'].append(j+1)
            output['specie'].append(key)
            output['wavenumber'].append(
                peaks[key]['R-Branch']['v_tilde(cm^-1)'][i]
            )
            output['abs'].insert(0, peaks[key]['R-Branch']['abs'][i])
            j += 1
        
        # transitions in p-branch go from 
        # highest wavenumber -> lowest wavenumber
        # fill P-Branch transitions
        j = 0
        for i in reversed(range(len(peaks[key]['P-Branch']['abs']))):
            peaks[key]['P-Branch']['j_values']['j_init'].insert(0, j+1)
            peaks[key]['P-Branch']['j_values']['j_fini'].insert(0, j)
            output['branch'].insert(0, 'P-Branch')
            output['j_init'].insert(0, j+1)
            output['j_fini'].insert(0, j)
            output['specie'].insert(0, key)
            output['wavenumber'].insert(
                0, 
                peaks[key]['P-Branch']['v_tilde(cm^-1)'][i]
            )
            output['abs'].insert(0, peaks[key]['P-Branch']['abs'][i])
            j += 1
    
    # sort by wavenumber 
    df = pd.DataFrame(output)
    df.sort_values(
        by=['wavenumber'], 
        ascending=True, 
        ignore_index=True,
        inplace=True)
    df.to_csv('transitions.csv')
    return df