#!/usr/bin/env python
# -*- coding:utf-8 -*-

#TODO: add dependence check

import os, subprocess

def build_clib():
    make_cmd = ['make']
    os.chdir(os.path.join(os.path.dirname(__file__), 'OneDQuantum'))
    print("Building C Lib")
    subprocess.check_call(['make'])
    try: 
        subprocess.check_call(['make', 'MP'])
    except CalledProcessError:
        print("openMP not supported")

if __name__ == "__main__":
    build_clib()

# vim: ts=4 sw=4 sts=4 expandtab
