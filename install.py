#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# TODO: add dependence check

import os, sys, subprocess


def build_clib(path, MSBuild=None):
    if not MSBuild:
        make_cmd = ['make']
        makemp_cmd = ['make', 'MP']
    else:
        make_cmd = [MSBuild, 'OneDQuantum.sln', '/p:Configuration=Release']
        makemp_cmd = [MSBuild, '1DSchrodinger.vcxproj',
                      '/p:Configuration=MP_Release']
        print(make_cmd)
    os.chdir(os.path.join(path, 'OneDQuantum'))
    print("Building C Lib")
    subprocess.check_call(make_cmd)
    try:
        subprocess.check_call(makemp_cmd)
    except CalledProcessError:
        print("openMP not supported")


def build_doc(path, MSBuild=None):
    if not MSBuild:
        make_cmd = ['make', 'html']
    else:
        print("MS building for documents is not now available")
        return
    os.chdir(os.path.join(path, 'docs'))
    print("Building Documents")
    try:
        subprocess.check_call(make_cmd)
    except subprocess.CalledProcessError:
        print("Building documents failed")


if __name__ == "__main__":
    MSBuild = None
    path = os.path.dirname(os.path.abspath(__file__))
    for opt in sys.argv[1:]:
        if opt.lower().startswith("msbuild="):
            MSBuild = opt[8:]
        else:
            print("Unknown option %s" % opt)
    build_clib(path, MSBuild)
    build_doc(path, MSBuild)

# vim: ts=4 sw=4 sts=4 expandtab
