#! /usr/bin/python3
#-*-coding: utf-8-*-

import numpy as np
import pandas as pd 

titanic = pd.read_csv("titanic.csv")
print(titanic)
print(titanic.head(8)) # to print the 8 first rows of the data frame
print(titanic.tail(10)) # to print the 10 last rows of the data frame
print(titanic.dtypes) # NO BRACKETS !!! dtypes is an attribute of a dataframe and series
#titanic.to_excel("titanic.xlsx", sheet_name = "passengers", index = False) # to_* methods used to store data => store the data as an excel file
# index = False => row index labels not saved in the spreadsheet
#titanic = pd.read_excel("titanic.xlsx", sheet_name = "passengers")
#print(titanic) # = error : install openpyxl

titanic.info() # to have a summary of the dataframe, technical info about a dataframe
# 891 entries = 891 rows
# Each row has a row label (aka the index) with values ranging from 0 to 890.
# The table has 12 columns. Most columns have a value for each of the rows (all 891 values are non-null). 
# Some columns do have missing values and less than 891 non-null values.
# The columns Name, Sex, Cabin and Embarked consists of textual data (strings, aka object). 
# The other columns are numerical data with some of them whole numbers (aka integer) and others are real numbers (aka float).
# The kind of data (characters, integers,…) in the different columns are summarized by listing the dtypes.
# The approximate amount of RAM used to hold the DataFrame is provided as well.

# REMEMBER
# Getting data in to pandas from many different file formats or data sources is supported by read_* functions.
# Exporting data out of pandas is provided by different to_*methods.
# The head/tail/info methods and the dtypes attribute are convenient for a first check.



###############################################
### Object creation
###############################################



s = pd.Series([1, 3, 5, np.nan, 6, 8])
#print(s)

dates = pd.date_range("20130101", periods = 6)
#print(dates)

df = pd.DataFrame(np.random.randn(6, 4), index = dates, columns = list("ABCD"))
#print(df)

df2 = pd.DataFrame(
    {
        "A" : 1.0,
        "B" : pd.Timestamp("20130102"),
        "C" : pd.Series(1, index = list(range(4)), dtype = "float32"),
        "D" : np.array([3] * 4, dtype="int32"),
        "E" : pd.Categorical(["test", "train", "test", "train"]),
        "F" : "foo",
    }
)
#print(df2)

#print(df2.dtypes)



###############################################
### Viewing data
###############################################



print(df.head())
print(df.tail(3))
print(df.index)
print(df.columns)


# WARNING !!! NumPy arrays have one dtype for the entire array, while pandas DataFrames have one dtype per column.
print(df.to_numpy()) 
print(repr(df.to_numpy())) # representation of the object : repr()
print(df2.to_numpy())
print(df.describe()) # to have a quick summary of the data
print(df.T) # to transpose the data
print(df.sort_index(axis = 1, ascending = False))
print(df.sort_values(by = "B"))



##############################################
### Getting
##############################################



print(df["A"])
print(df[0:3])
print(df["20130102":"20130104"])



##############################################
### Selection by label
##############################################




print(df.loc[dates[0]]) # Pandas provide a unique method to retrieve rows from a Data frame. 
DataFrame.loc[] method is a method that takes only index labels and returns row or dataframe if the index label exists in the caller data frame.

print(df.loc[:, ["A", "B"]])
print(df.loc["20130102":"20130104", ["A", "B"]])
print(df.loc["20130102", ["A", "B"]])
print(df.loc[dates[0], "A"]) # scalar value : index 0 of the dates, column "A"



############################################
### Selection by position
############################################