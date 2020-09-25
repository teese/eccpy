from time import strftime
from typing import Union

import pytest
import eccpy
from tempfile import TemporaryDirectory
from pathlib import Path

from test.helpers.helpers import set_up_excel_settings_for_functional_test, get_replacement_dict_with_test_paths


@pytest.mark.regression
def test_processing_of_example_data():
    with TemporaryDirectory() as tmpdirname:
        # uncomment to view test output in repo dir during development
        # tmpdirname: Union[Path, str] = eccpy_module_path / "test"

        tmpdir: Path = Path(tmpdirname)
        temp_path_for_testing: Union[Path, str] = tmpdir / "eccpy_temp_test_output"
        print(f" Test output will be saved to {temp_path_for_testing}, and deleted after tests are finished")

        replacement_dict: dict = get_replacement_dict_with_test_paths(temp_path_for_testing)
        test_settings: Union[Path, str] = set_up_excel_settings_for_functional_test(temp_path_for_testing, replacement_dict)

        eccpy.run_curvefit(test_settings)
        eccpy.run_gatherer(test_settings)

        assert_curvefit_output_files_exist(temp_path_for_testing)

        assert_gatherer_output_files_exist(temp_path_for_testing)

    # assert temporary output files are removed
    assert not Path(temp_path_for_testing).is_dir()


def assert_gatherer_output_files_exist(temp_path_for_testing):
    date_string = strftime("%Y%m%d")
    assert (temp_path_for_testing / f"analysed/{date_string}").is_dir()
    assert (temp_path_for_testing / f"analysed/{date_string}/{date_string}.xlsx").is_file()
    assert (temp_path_for_testing / f"analysed/{date_string}/{date_string}_datapoints.png").is_file()


def assert_curvefit_output_files_exist(temp_path_for_testing):
    assert (temp_path_for_testing / "analysed").is_dir()
    assert (temp_path_for_testing / f"generated_data_0/curves/AA control.png").is_file()
    assert (temp_path_for_testing / f"generated_data_0/generated_data_0_summary.png").is_file()
