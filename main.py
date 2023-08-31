'''
Created on Oct 1, 2020

@author: Ryan
'''
import matplotlib
matplotlib.interactive(True)  # show non-blocking
matplotlib.use('WXAgg')
import wx
import pathlib
import argparse

from gcms_align import main_ui


INTERACTIVE_MULTIPROCESS = True
VERSION = "1.0"


if __name__ == '__main__':
    print("STATUS: Minerva Starting Setup")
    app = wx.App(0)
    if not INTERACTIVE_MULTIPROCESS:
        parser = argparse.ArgumentParser(description='Analyze GCMS files')
        parser.add_argument("-d", "--target_directory", help="Directory containing data")
        parser.add_argument("-c", "--config_file", help="Text file for configuration settings")
        parser.add_argument("-p", "--parse_only", help="Disable user interface", action="store_true")
        parser.add_argument("--version", action="version", version=f"{VERSION}")
        args = parser.parse_args()
        if args.target_directory:
            input_dir = pathlib.Path(args.target_directory)
        else:
            input_dir = None
        if args.config_file:
            config_file = pathlib.Path(args.config_file)
        else:
            config_file = None

        if args.parse_only:
            main_ui.MainProcessing(version=VERSION, config_file=config_file, data_input_dir=input_dir)
    else:
        ui = main_ui.MainUI(version=VERSION, default_input_dir=None)
        print("STATUS: Setup complete, displaying UI")
        app.MainLoop()
