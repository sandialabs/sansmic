# coding: utf-8
#
# SPDX-License-Identifier: MIT-License
# 
# Copyright (c) 2020 Aubrey Taylor
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

def print_progress(iteration, total, prefix="", suffix="", decimals=1, length=100, fill="█", printEnd="\r", useRatio=False):
    r"""Call in a loop to create terminal progress bar.

    Source: https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a

    Parameters
    ----------
    iteration : int
        current iteration
    total : int
        total iterations
    prefix : str, optional
        prefix string for the start of the progress line, by default ''
    suffix : str, optional
        suffix string to go at the end of the progress line, by default ''
    decimals : int, optional
        positive number of decimal places to include in the percentage count, by default 1
    length : int, optional
        the number of characters in the progress bar, by default 100
    fill : str, optional
        the character to use for completed progress, by default ``█``
    printEnd : str, optional
        the character to use at the end of the line, by default ``'\r'`` (carriage return without newline)
    useRatio : bool, optional
        whether the progress should be printed as a percentage (the default) or as a ratio, by default False
    """
    nc = len(str(total))
    ratio = ("{:>" + str(nc) + "d}/{:<" + str(nc) + "d}").format(iteration, total)
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    if not useRatio:
        print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    else:
        print(f"\r{prefix} |{bar}| {ratio} {suffix}", end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()
