from getpass import getuser as uid
from datetime import datetime as dt
from custom_types import *

def box(string: str, delim_char: char) -> str:
    outline = len(string) + 4
    spacer = delim_char * outline
    top = f'{spacer:<{outline}}'
    middle = (f'\n{f'{delim_char} {string} {delim_char}':<{outline}}')
    bottom = (f'\n{spacer:^{outline}}')
    box = top + middle + bottom 
    return box

class outputData():
    """
    Tracks data from analysis 
    """
    
    def __init__(self, save_path: str, specie: str):
        
        """
        Initializes output. Run at beginning of script execution.
        """
        
        bbox = f'{box(specie, '@')}'
        user = f'USER:\t{uid()}'
        date = f'DATE:\t{dt.now().date()}'
        time = f'TIME:\t{dt.now().time()} ({dt.now().astimezone().tzinfo})'
        expr = f'EXPR:\t{specie}'
        _id_ = f'{user}\n{date}\n{time}\n{expr}\n'
        self.__lines = [_id_, bbox]
        self.__cwd = save_path
        self.__specie = specie
        
    def __add_lines(self, line):
        self.__lines.append(line)
    
    def StateVars(self, decimals: int, state_vars: dict):
        
        specie = self.__specie
        
        bbox = box(f'State Variables for {specie}', '*')
        self.__lines.append(bbox)
        
        fmt = f".{decimals}f"
        
        # left column (state variable name)
        l_col_width = 32
        
        # center column (value)
        c_col_width = 7 + decimals        
        
        # right column (units)
        r_col_width = 6
        
        # column titles
        l_col = f'{f'State Variable':<{l_col_width}}'
        c_col = f'{f'Value':<{c_col_width}}'
        r_col = f'{f'Units':<{r_col_width}}'
        spacer = '-' * ( 
                        l_col_width
                        + c_col_width
                        + r_col_width
                        )
        out = f'{spacer}\n{l_col}{c_col}{r_col}\n'
        for var in state_vars.keys():
            
            # skip any kg masses (too small for rounding)
            if var == 'spring-const':
                var_name = 'Spring Constant (k)'
            elif var == 'mass1-amu':
                var_name = 'Mass 1'
            elif var == 'mass2-amu':
                var_name = 'Mass 2'
            elif var == 'reduced-mass-amu':
                var_name = 'Reduced Mass (μ)'
            else:
                continue
            
            value = format(state_vars[var][0], fmt)
            unit = state_vars[var][1]
            state_var_str = f'{f'{var_name}':<{l_col_width}}'
            value_str = f'{f'{value}':<{c_col_width}}'
            unit_str = f'{f'{unit}':<{r_col_width}}'
            out += state_var_str + value_str + unit_str
            out += '\n'
        
        self.__lines.append(out)
        
    def Const(
            self,
            analysis_name: str,
            titles: nd_1DStrArray,
            decimals: int,
            r_v: nd_1DArray,
            p_v: nd_1DArray,
            avg: nd_1DArray, 
            se: nd_1DArray,
            units: nd_1DStrArray
        ):
        
        
        # left column (constant name)
        l_col_width = 32
        
        # left middle column (P-Branch)
        lc_col_width = 7 + decimals
        
        # center middle column (R-Branch)
        mc_col_width = 7 + decimals        
        
        # right middle column (avg value)
        rc_col_width = 16 + decimals   
        
        # right column (units)
        r_col_width = 6
        
        # column titles
        l_col = f'{f'Constant':<{l_col_width}}'
        lc_col = f'{f'P-Branch':<{lc_col_width}}'
        mc_col = f'{f'R-Branch':<{mc_col_width}}'
        rc_col = f'{f'Average':<{rc_col_width}}'
        r_col = f'{f'Units':<{r_col_width}}'
        spacer = '-' * (
                        l_col_width 
                        + lc_col_width
                        + mc_col_width
                        + rc_col_width 
                        + r_col_width
                        )
        out = f'{spacer}\n{l_col}{lc_col}{mc_col}{rc_col}{r_col}\n'
        
        for i in range(len(titles)):
            # left column: titles
            title_str = f'{titles[i]}:'
            title_out =  f'{title_str:<{l_col_width}}'
            out += title_out
            
            fmt = f".{decimals}f"
            
            # left center column: p-branch values
            p_val_str = f'{format(p_v[i], fmt)}'
            p_val_out = f'{p_val_str:<{lc_col_width}}'
            out += p_val_out
            
            # middle center column: r-branch values
            r_val_str = f'{format(r_v[i], fmt)}'
            r_val_out = f'{r_val_str:<{mc_col_width}}'
            out += r_val_out
            
            # right center column: average values
            avg_fmt = f'{format(avg[i], fmt)}'
            se_fmt = f'{format(se[i], fmt)}'
            value_str = f'({avg_fmt}±{se_fmt})'
            value_out = f'{value_str:<{rc_col_width}}'
            out += value_out
            
            # right column: units
            unit_out =  f'{units[i]:<{r_col_width}}'
            out += unit_out
            out += '\n'
        
        bbox = box(analysis_name, '*')
        self.__add_lines(f'{bbox}\n{out}')
    
    def output(self):
        f_name = f'{self.__cwd}/{self.__specie.lower()}.txt'
        with open(f_name, 'w') as f:
            for line in self.__lines:
                f.write(f'{line}\n\n')