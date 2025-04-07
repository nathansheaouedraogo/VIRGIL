import pandas as pd

def analyze_branch(specie: str, branch_name: str) -> tuple[float, float]:
    while True:
        branch_start = float(input(
                    f'\nPlease input the start of the {branch_name} for {specie}: '
                    ))
        branch_end = float(input(
                    f'\nPlease input the end of the {branch_name} for {specie}: '
                    ))
        out_1 = f'\n{branch_name} range for {specie}: '
        out_2 = f'\t{branch_start}cm^(-1)  to {branch_end}cm^(-1)'
        print(out_1+out_2)
        prompt = input(f'\nis this correct? [y/n]: ')
        if prompt == 'y':
            print('\n')
            return branch_start, branch_end
        elif prompt == 'n':
            continue
        else:
            print('\nexiting...')
            exit(1)


def get_bounds(branch_name: str) -> tuple[float, float]:
    while True:    
        inst_1 = f'\nPlease set minimum abs threshold for peaks: '
        threshold = float(input(inst_1))
        
        inst_1 = f'Please set a minimum spacing between peaks: '
        spacing = float(input(inst_1))

        out_1 = f'\nMinimum peak threshold for {branch_name}: {threshold} abs'
        out_2 = f'\nMinimum peak spacing   for {branch_name}: {spacing} cm^(-1)'

        print(out_1+out_2)
        prompt = input(f'\nis this correct? [y/n]: ')
        if prompt == 'y':
            print('\n')
            return threshold, spacing
        elif prompt == 'n':
            continue
        else:
            print('\nexiting...')
            exit(1)


def find_peaks(
    df_data: pd.DataFrame, 
    x_start: float, 
    x_end: float, 
    y_threshold: float, 
    x_name: str, 
    y_name: str, 
    x_spacing: float = 0
    ) -> pd.DataFrame:
    
    """
    Finds values in y which are above minimum threshold in a dataframe. 
    >>> NOTE: 
            In order to properly handle edge cases (first, last peaks),
            x_start should be slightly before the first peak 
            and x_end should be slightly after the second peak
    Args:
        df_data (pd.DataFrame): dataframe containing x and y values 
        x_start (float): start of data range
        x_end (float): end of data range
        y_threshold (float): any values >= to this value will be processed
        x_name (str): name of x column in df_data
        y_name (str): name of y column in df_data
        x_spacing (float, optional): minimum spacing between peaks. Defaults to 0.

    Returns:
        pd.DataFrame: maxima (from highest x values -> lowest x value)
    """
            
    # 1. filter branch
    
    mask = (
        (df_data[x_name] >= x_start) &
        (df_data[x_name] <= x_end) &
        (df_data[y_name] >= y_threshold)
    )
    
    df_maxima = df_data[mask].copy()
    
    # 2. Handle special cases: empty DataFrame, negative x_spacing, or no spacing filter.
    if df_maxima.empty:
        output = df_maxima
    elif x_spacing < 0:
        print(f'\nWarning: x_spacing must be >= 0 (got: {x_spacing})... ignoring spacing')
        output = df_maxima
    elif x_spacing == 0:
        output = df_maxima
    else:
        # 3. Convert candidate peak columns to lists for easy iteration.
        y_maxima = df_maxima[y_name].to_list()  # List of y-values (amplitudes).
        x_maxima = df_maxima[x_name].to_list()  # List of x-values.
        i_maxima = df_maxima.index.to_list()    # Original indices of candidate peaks.
        
        # 4. Initialize the final accepted indices with the first candidate.
        # ensures we keep the first peak in the order of appearance.
        f_maxima = [i_maxima[0]]
        last_accepted_x = x_maxima[0]
        last_accepted_y = y_maxima[0]
        
        # 5. Iterate over remaining candidate peaks.
        for j in range(1, len(i_maxima)):
            current_x = x_maxima[j]
            current_y = y_maxima[j]
            
            # check if the current peak is too close to the last accepted peak.
            if abs(current_x - last_accepted_x) < x_spacing:
                
                # if the current peak's amplitude is higher, replace the last accepted peak.
                if current_y > last_accepted_y:
                    f_maxima[-1] = i_maxima[j]
                    last_accepted_x = current_x
                    last_accepted_y = current_y
                
                # Else, retain the previously accepted peak (do nothing).
            else:
                # If the current peak is sufficiently distant, accept it.
                f_maxima.append(i_maxima[j])
                last_accepted_x = current_x
                last_accepted_y = current_y
                
        # 6 Filter the candidate DataFrame to include only the accepted peaks.
        mask = ~df_maxima.index.isin(f_maxima)
        df_maxima.drop(df_maxima[mask].index, inplace=True)
        df_maxima.reset_index(inplace=True, drop=True)
        
        # NOTE: 7 -> 9 deprecated; usr must now ensure x_start < first peak and x_end > last peak
        # # 7. check if peaks @ x_start meets threshold
        # #    if so, process 
        # init_x = x_start
        # init_y = df_data[y_name][df_data[x_name] == init_x].values[0]
        # process_init = (init_y != df_maxima[y_name].iloc[-1]) and (init_y >= y_threshold)
        # if process_init:
        #     curr_init_x = df_maxima[x_name].iloc[-1]
        #     curr_init_y = df_maxima[y_name].iloc[-1]
            
        #     # if dst < x_spacing, 
        #     # change value at idx[-1] iff y1 > df[y col].idx[-1]
        #     if abs(init_x-curr_init_x) < x_spacing:
        #         if init_y > curr_init_y:
        #             df_maxima[x_name].iat[-1] = init_x
        #             df_maxima[y_name].iat[-1] = init_y
            
        #     # else, prepend to df and reset axis
        #     else:
        #         init_line = pd.DataFrame({
        #             f'{x_name}' : init_x,
        #             f'{y_name}' : init_y
        #         }, index=[0])
        #         df_maxima = pd.concat(
        #             [init_line, df_maxima.iloc[:,:]]
        #         ).reset_index(drop=True)
        
        
        # # 8. check if peaks @ x_start meets threshold
        # #    if so, process 
        # fini_x = x_end
        # fini_y = df_data[y_name][df_data[x_name] == fini_x].values[0]
        # process_fini = (fini_x != df_maxima[y_name].iloc[0]) and (fini_y >= y_threshold)
        # if process_fini:
        #     curr_fini_x = df_maxima[x_name].iloc[0]
        #     curr_fini_y = df_maxima[y_name].iloc[0]
            
        #     # if dst < x_spacing, 
        #     # change value at idx[0] iff y1 > df[y col].idx[0]
        #     if abs(fini_x-curr_fini_x) < x_spacing:
        #         if fini_y > curr_fini_y:
        #             df_maxima[x_name].iat[0] = fini_x
        #             df_maxima[y_name].iat[0] = fini_y
            
        #     # else, prepend to df and reset axis
        #     else:
        #         fini_line = pd.DataFrame({
        #             f'{x_name}' : init_x,
        #             f'{y_name}' : init_y
        #         }, index=[0])
        #         df_maxima = pd.concat(
        #             [df_maxima.iloc[:,:], fini_line]
        #         ).reset_index(drop=True)
                
        # 9. since peaks are reversed (lowest -> highest), 
        # inverse dataframe before outputting
        
        # convert to numeric values
        df_maxima[x_name] = pd.to_numeric(df_maxima[x_name])
        df_maxima[y_name] = pd.to_numeric(df_maxima[y_name])
        output = df_maxima[::-1].reset_index(drop=True)

    # 10. output peaks
    return output 
        
def peak_pick(
    peaks : dict, 
    df : pd.DataFrame,
    show_plot: bool,
    x_name: str, 
    y_name: str
    ) -> None:
    
    if show_plot == True:
        from interactive_plot import plot_peaks
        plot_peaks(df[x_name], df[y_name])

    
    for key in peaks.keys():

        # analyze p-branch
        p_start, p_end, = analyze_branch(peaks[key]['specie'], 'P-Branch')
        thresh, spacing = get_bounds('P-Branch')
        df_pBranch = find_peaks(
            df, p_start, p_end, thresh, x_name, y_name, spacing
            )
        peaks[key]['P-Branch']['abs'].extend(
            df_pBranch['abs'].to_list())
        peaks[key]['P-Branch']['v_tilde(cm^-1)'].extend(
            df_pBranch['v_tilde(cm^-1)'].to_list())
        
        # analyze p-branch
        r_start, r_end, = analyze_branch(peaks[key]['specie'], 'R-Branch')
        thresh, spacing = get_bounds('R-Branch')
        df_rBranch = find_peaks(
            df, r_start, r_end, thresh, x_name, y_name, spacing
            )
        peaks[key]['R-Branch']['abs'].extend(
            df_rBranch['abs'].to_list())
        peaks[key]['R-Branch']['v_tilde(cm^-1)'].extend(
            df_rBranch['v_tilde(cm^-1)'].to_list())

        # print newline buffer
        print("\n\n")