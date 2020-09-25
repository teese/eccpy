import pytest
from eccpy.tools import assert_df_contains_no_nan_values
import pandas as pd
import numpy as np

@pytest.mark.regression
def test_assert_df_contains_no_nan_values():

    data = [[np.nan,2],[3,np.nan]]
    df = pd.DataFrame(data)
    df.columns = ["col1", "col2"]
    df.index = ["ind1", "ind2"]

    with pytest.raises(ValueError) as ve:
        assert_df_contains_no_nan_values(df)

    data = [[1,2],[3,4]]
    df = pd.DataFrame(data)
    assert_df_contains_no_nan_values(df)

