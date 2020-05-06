#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# TODO: add dependence check

import os, sys, subprocess
from subprocess import CalledProcessError


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
    if MSBuild is None:
        make_cmd = ['make', 'html']
    else:
        print("MS building for documents is not now available")
        return
    #  os.chdir(os.path.join(path, 'OneDQuantum/docs'))
    os.chdir(os.path.join(path, 'docs'))
    print("Building Documents")
    try:
        subprocess.check_call(make_cmd)
    except subprocess.CalledProcessError:
        print("Building documents failed")


if __name__ == "__main__":
    MSBuild = None
    path = os.path.dirname(os.path.abspath(__file__))
    docs = None
    for opt in sys.argv[1:]:
        if opt.lower().startswith("--msbuild="):
            MSBuild = opt[10:]
        elif opt.lower().startswith("--nodoc"):
            docs = False
        elif opt.lower().startswith("--doc"):
            docs = True
        else:
            print("Unknown option %s" % opt)
    build_clib(path, MSBuild)
    if docs is None:
        if input("Build the documentation locally? Y/[N]\n"
                 ).lower().startswith('y'):
            docs = True
    if docs:
        build_doc(path, MSBuild)

# vim: ts=4 sw=4 sts=4 expandtab
