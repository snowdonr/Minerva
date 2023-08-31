'''
Created on Nov. 14, 2021

@author: Ryan
'''
import subprocess
import pathlib
import shutil
import os


if __name__ == '__main__':
    pyinstall_path = r"C:\Python37\Scripts\pyinstaller.exe"
    minerva_path = r"D:\Ryan\University of Calgary\STRIDE - Minerva - Minerva\Source\src\main.py"
    icon_info = r'--icon=D:\Ryan\University of Calgary\STRIDE - Minerva - Minerva\Minerva_Book.ico'
    hidden_import = r"--hidden-import"
    module_cf = "cftime"
    module_pymz = "pymzml"
    confirm_overwrite = "-y"
    output_dir = r"C:\Build_Minerva"
    subprocess.run(
        [pyinstall_path, minerva_path, icon_info, hidden_import, module_cf, hidden_import, module_pymz, confirm_overwrite],
        cwd=output_dir,
        check=True  # , capture_output=True
    )

    # Copy pymzml and pyms_nist_search from C:\Python37\Lib\site-packages to dist\main directory
    package_source = pathlib.Path(r"C:\Python37\Lib\site-packages")
    package_names = ["pymzml", "pyms_nist_search/templates"]
    target_dir = pathlib.Path(output_dir)/"dist"/"main"

    for package in package_names:
        shutil.copytree(package_source/package, target_dir/package)  # , dirs_exist_ok=True)
    #
    # rename main to Minerva
    os.rename(target_dir, target_dir.parent/"Minerva")

    # TODO: Create zipfile?
    # TODO: check dir, create inno config
    # TODO: run inno

    print("Finished")
    # Ubuntu packages
    # scipy, pandas, matplotlib, wxPython, libgtkmm-3.0-1, libnotify4, sqlalchemy, netcdf4, pymzml, openpyxl
    # pyms-nist-search
    # biopython  TODO used only for treecluster
