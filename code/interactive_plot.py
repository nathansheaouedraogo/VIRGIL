import pandas as pd
from custom_types import *
import plotly.express as px

def plot_peaks(
    x: pd.Series, 
    y: pd.Series
    ) -> px.Figure:
    
    plot = px.line(
        x = x, 
        y =y,
        labels = {
            'x' : 'wavenumber (cm^-1)',
            'y' : 'Absorbance'
        }
    )
    return plot