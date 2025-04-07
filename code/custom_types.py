from typing import Annotated, Literal
import numpy as np
import numpy.typing as npt

# 1D array of floats: expected shape (n,)
nd_1DArray: Annotated[
    npt.NDArray[np.float64], 
    Literal[1]
    ] = np.array([1.0, 2.0, 3.0])

# 1D array of strings
nd_1DStrArray: Annotated[
    npt.NDArray[np.string_], 
    Literal[1]
    ] = np.array(['str1', 'str2', 'str3'])


# Character type that ensures only one character
class char(str):
    def __new__(cls, s: str) -> "char":
        if len(s) != 1:
            raise ValueError("Only one character allowed")
        return super().__new__(cls, s)

