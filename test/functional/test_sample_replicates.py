import glob
import pytest
import eccpy
import pandas as pd
from eccpy.tools import get_eccpy_module_path
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union
from shutil import copyfile

from test.helpers.helpers import get_replacement_dict_with_test_paths, set_up_excel_settings_for_functional_test, replace_xlsx_cells_by_replacement_dict


@pytest.mark.regression
def test_processing_of_example_data():
    eccpy_module_path = get_eccpy_module_path()

    with TemporaryDirectory() as tmpdirname:
        # uncomment to view test output in repo dir during development
        # tmpdirname: Union[Path, str] = eccpy_module_path / "test"

        tmpdir: Path = Path(tmpdirname)
        temp_path_for_testing: Union[Path, str] = tmpdir / "eccpy_temp_test_output"
        print(f" Test output will be saved to {temp_path_for_testing}, and deleted after tests are finished")

        example_data_dir: Path = eccpy_module_path / "eccpy/examples/example_data"
        data_files_orig: list = glob.glob(str(example_data_dir / "*.xlsx"))
        temp_example_data_dir: Path = temp_path_for_testing / "examples"
        for data_file_path in data_files_orig:
            assert Path(data_file_path).is_file()
            destination = temp_example_data_dir / Path(data_file_path).name
            if not destination.parent.is_dir():
                destination.parent.mkdir(parents=True)
            copyfile(data_file_path, destination)

        data_file_0 = temp_example_data_dir / "generated_data_0.xlsx"
        dose_replacements = {"B1": "sample1", "C1": "sample1", "D1": "sample1"}
        response_replacements = {"B1": "sample1", "C1": "sample1", "D1": "sample1"}

        replacement_dict = {}
        replacement_dict["dose"] = dose_replacements
        replacement_dict["response"] = response_replacements

        replace_xlsx_cells_by_replacement_dict(replacement_dict, data_file_0)

        replacement_dict: dict = get_replacement_dict_with_test_paths(temp_path_for_testing, temp_example_data_dir)
        test_settings: Union[Path, str] = set_up_excel_settings_for_functional_test(temp_path_for_testing, replacement_dict)

        eccpy.run_curvefit(test_settings)
        eccpy.run_gatherer(test_settings)

        output_data_0_xlsx: Path = temp_path_for_testing / f"generated_data_0/generated_data_0.xlsx"
        assert output_data_0_xlsx.is_file()

        df = pd.read_excel(output_data_0_xlsx, sheet_name="by_sample", index_col=0)
        total_n_sample1_in_by_sample_tab = df.at["sample1", "n"] + df.at["sample1", "n_excluded"]
        assert total_n_sample1_in_by_sample_tab == len(dose_replacements)
        assert "mean" in df.columns
        assert "std" in df.columns

    # assert temporary output files are removed
    assert not Path(temp_path_for_testing).is_dir()
