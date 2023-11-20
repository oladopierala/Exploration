#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 11:35:37 2023

@author: giuliaserino
"""


import os
import glob
import pandas as pd

# Define the folder path where the videocoding files (txt) are located
folder_path = r'/Users/giuliaserino/Library/CloudStorage/OneDrive-Birkbeck,UniversityofLondon/Canada/video_coding/export//'
# Define the output folder path where the videocoding files (xlsx) will be located
output_path = r'/Users/giuliaserino/Library/CloudStorage/OneDrive-Birkbeck,UniversityofLondon/Canada/video_coding/SA//'

# Define the column names
column_names = ['event', 'n','start_time', 'end_time','duration', 'l']

# Get a list of all text files in the folder
file_list = glob.glob(os.path.join(folder_path, '*.txt'))

# Iterate over each file
for file_path in file_list:
    # Extract the file name without the extension
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Read the DataFrame from the file
    df = pd.read_csv(file_path, header=None, delimiter='\t')
    
    
    # Set the column names
    df.columns = column_names
    
    # select column of interest 
    df = df[['event', 'start_time', 'end_time','duration']]
    
    # Define the output CSV file path using the same name as the text file
    output_csv_path = os.path.join(output_path, file_name + '.xlsx')

    # Export the DataFrame to a xlsx file
    
    df.to_excel(output_csv_path, index=False)