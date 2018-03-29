#!/usr/bin/python
import re


def check_delimiter(sbtab_file):
    '''
    determine the delimiter of the tabular file
    '''
    sep = False

    for row in sbtab_file.split('\n'):
        if row.startswith('!!'): continue
        if row.startswith('!'):
            s = re.search('(.)(!)', row[1:])
            # if there is only 1 column, we have to define a default separator
            # let's use a tab.
            try: sep = s.group(1)
            except: sep = '\t'

    return sep
