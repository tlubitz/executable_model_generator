#!/usr/bin/python
import re
import xlwt


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


def csv2xls(sbtabs, sbtab_types, filename):
    '''
    converts sbtab file to xls file
    @sbtab_file: sbtab string
    '''
    book = xlwt.Workbook()

    if len(sbtabs) > 1:
        for k, sbtab in enumerate(sbtabs):
            sheet = book.add_sheet(sbtab_types[k].capitalize())

            header_row = sheet.row(0)
            header_row.write(0, sbtab.header_row)

            column_row = sheet.row(1)
            for i, column in enumerate(sbtab.columns):
                column_row.write(i, column)

            for i, row in enumerate(sbtab.value_rows):
                xls_row = sheet.row(i + 2)
                for j, entry in enumerate(row):
                    xls_row.write(j, entry)

    if filename.endswith('.xls'): book.save(filename)
    else: book.save(filename + '.xls')   
