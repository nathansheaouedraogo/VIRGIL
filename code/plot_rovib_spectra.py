import pandas as pd
import matplotlib.pyplot as plt

def plot_rovib_spectra(
    df: pd.DataFrame, 
    peaks: dict, 
    show_fig: bool  
    ) -> plt.Figure:
    
    # plot!
    fig, axis = plt.subplots()
    axis.plot(
        df['v_tilde(cm^-1)'], 
        df['abs'],
        color = 'black',
        linewidth = 1
    )
    axis.set_xlim(left=df['v_tilde(cm^-1)'].min(), right=df['v_tilde(cm^-1)'].max()) 
    axis.set_xlabel(
        rf'$\tilde{{\nu}} \quad \left(\mathrm{{cm}}^{{-1}}\right)$',
        fontsize=16
        )
    axis2 = axis.secondary_xaxis('top')
    x_ticks = [
        peaks[key]['constants']['nuTildeNaught']['estimated'][0][0] for key in peaks.keys()
        ]
    axis2.set_xticks(x_ticks)
    axis2.set_xticklabels([f'{round(x)}' for x in x_ticks])

    axis2.set_xlabel(
        rf'$\text{{Estimated}}\quad \tilde{{\nu}}_{{0}} \quad \text{{Bands}} \quad \left(\mathrm{{cm}}^{{-1}}\right)$',
        fontsize=14
    )

    axis.set_ylabel('Absorbance', fontsize=14)
    y_max = 0.94
    axis.set_ylim(bottom=0, top=1)


    # branches and fundamental vib band
    for key in peaks.keys():  
        
        
            
        # P-Branch
        p_branch = axis.axvspan(
            xmin = peaks[key]['P-Branch']['v_tilde(cm^-1)'][0],
            xmax = peaks[key]['P-Branch']['v_tilde(cm^-1)'][-1],
            alpha = 0.25,
            color = peaks[key]['colour']['branch']
        )
        axis.text(
            x = (
                peaks[key]['P-Branch']['v_tilde(cm^-1)'][0]
                + peaks[key]['P-Branch']['v_tilde(cm^-1)'][-1]
                ) / 2, 
            y = y_max,
            s = "P-Branch",
            ha = 'center', 
            va = 'center', 
            fontsize = 12,
        )
        
        # Q-Branch
        q_branch = axis.axvspan(
            xmin = peaks[key]['P-Branch']['v_tilde(cm^-1)'][-1],
            xmax = peaks[key]['R-Branch']['v_tilde(cm^-1)'][0],
            alpha = 0.5,
            color = peaks[key]['colour']['center']
        )
        # # vib band
        # axis.axvline(
        #     ymax = y_max,
        #     linestyle = 'dashdot',
        #     color='black',
        #     linewidth = 2
        #     )
        axis.text(
            x = (
                peaks[key]['P-Branch']['v_tilde(cm^-1)'][-1]
                + peaks[key]['R-Branch']['v_tilde(cm^-1)'][0]
                ) / 2, 
            y = y_max,
            rotation = 90,
            s = "Q-Branch",
            ha = 'center', 
            va = 'top', 
            fontsize = 8,
        )
        
        # R-Branch
        r_branch = axis.axvspan(
            xmin = peaks[key]['R-Branch']['v_tilde(cm^-1)'][0],
            xmax = peaks[key]['R-Branch']['v_tilde(cm^-1)'][-1],
            alpha = 0.25,
            color = peaks[key]['colour']['branch']
        )
        axis.text(
            x = (
                peaks[key]['R-Branch']['v_tilde(cm^-1)'][0]
                + peaks[key]['R-Branch']['v_tilde(cm^-1)'][-1]
                ) / 2, 
            y = y_max,
            s = "R-Branch",
            ha = 'center', 
            va = 'center', 
            fontsize = 12,
        )
        
        # label info, vib band
        l_1 = rf'$\tilde{{\nu}}_{{0,\text{{{peaks[key]['specie']}}}}}$'
        l_2 = rf'$\approx {round(peaks[key]['constants']['nuTildeNaught']['estimated'][0][0])}\ \mathrm{{cm}}^{{-1}}$' 
        label = l_1 + l_2 
        
        axis.text(
                x = peaks[key]['constants']['nuTildeNaught']['estimated'][0][0], 
                y = -0.1, 
                s = label,
                ha='center', 
                va='bottom', 
                fontsize = 12, 
                bbox=dict(
                    facecolor = peaks[key]['colour']['center'], 
                    alpha = 0.5, 
                    boxstyle = 'round,pad=0.5'
                )
            )
    
    fig.tight_layout()
    
    if show_fig == True:
        plt.show()
    
    return fig 

