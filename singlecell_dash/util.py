import numpy as np
import pandas as pd


def maybe_format(item):
    """Pretty-format a string, integer, float, or percent

    Parameters
    ----------
    item : pandas.Series
        A single-item series containing a .name attribute and a value in the
        first (0th) index
    """
    value = item[0]
    if pd.isnull(value):
        return 'N/A'
    elif isinstance(value, str):
        return value
    elif 'percent' in item.name.lower():
        return '{:.2f}%'.format(value)
    elif isinstance(value, pd.Timestamp):
        return str(np.datetime64(value, 'D'))
    elif (isinstance(value, float)  # this must go before ints!
          or np.issubdtype(value, np.number)):
        if value >= 1e3:
            return locale.format("%d", int(value), grouping=True)
        else:
            return locale.format("%.3g", value, grouping=True)
    elif (isinstance(value, int)
          or np.issubdtype(value, np.integer)):
        return locale.format("%d", value, grouping=True)
    else:
        raise TypeError