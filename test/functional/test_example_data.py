from typing import Union

import pytest
import eccpy
from eccpy.tools import get_eccpy_module_path
from openpyxl import load_workbook
from pathlib import Path
from shutil import copyfile, rmtree


@pytest.mark.regression
def test_processing_of_example_data():
    eccpy_module_path: Union[Path, str] = get_eccpy_module_path()
    test_temp_output_path: Union[Path, str] = eccpy_module_path / "test/temp_output"

    test_settings: Union[Path, str] = set_up_excel_settings_for_functional_test()

    eccpy.run_curvefit(test_settings)
    eccpy.run_gatherer(test_settings)

    rmtree(test_temp_output_path)



def set_up_excel_settings_for_functional_test():
    eccpy_module_path = get_eccpy_module_path()

    orig_settings_xlsx: Union[Path, str] = eccpy_module_path / "docs/example_settings/ECCpy_settings_template.xlsx"
    temp_test_settings_xlsx: Union[Path, str] = eccpy_module_path / "test/temp_output/settings_example_data.xlsx"

    example_data_dir: str = str(eccpy_module_path / "docs/example_data")
    output_data_dir: str = str(eccpy_module_path / "test/temp_output")

    if not temp_test_settings_xlsx.parent.is_dir():
        temp_test_settings_xlsx.parent.mkdir()

    copyfile(orig_settings_xlsx, temp_test_settings_xlsx)
    assert temp_test_settings_xlsx.is_file()

    wb = load_workbook(temp_test_settings_xlsx)
    sheet = "files"
    input_data_cells = ["G2", "G3", "G4", "G5"]
    output_data_cells = ["H2", "H3", "H4", "H5"]

    ws = wb.get_sheet_by_name(sheet)

    for cell in input_data_cells:
        ws[cell] = example_data_dir
    for cell in output_data_cells:
        ws[cell] = output_data_dir

    wb.save(temp_test_settings_xlsx)

    return temp_test_settings_xlsx


