'''
Created on May 6, 2022

@author: Ryan
'''
from gcms import cirpy
# import openpyxl
import logging
import pandas as pd


class LibraryUtil(object):

    def __init__(self, excel_path):
        if not excel_path.is_file():
            logging.warning(f"Could not open excel file {excel_path}")
            return

        df_compound = pd.DataFrame(
            pd.read_excel(excel_path, sheet_name="Compounds")
        )
        
        # wb = openpyxl.load_workbook(excel_path)

        result_list = cirpy.resolve("ethanol".casefold(), "cas")
        print(result_list)


if __name__ == '__main__':
    base = LibraryUtil()
