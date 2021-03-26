#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# TODO: add dependence check
import argparse

import os
import sys
import subprocess
from subprocess import CalledProcessError
from create_shortcut import create_shortcut


def build_clib(path, MSBuild=None):
    if MSBuild is None:
        make_cmd = ['make']
        makemp_cmd = ['make', 'MP']
    else:
        make_cmd = [MSBuild, 'OneDQuantum.sln', '/p:Configuration=Release']
        makemp_cmd = [MSBuild, '1DSchrodinger.vcxproj',
                      '/p:Configuration=MP_Release']
    os.chdir(os.path.join(path, 'OneDQuantum'))
    print("Building C Lib")
    subprocess.check_call(make_cmd)
    try:
        subprocess.check_call(makemp_cmd)
    except CalledProcessError:
        print("openMP not supported")


def build_doc(path):
    if sys.platform.startswith('win'):
        print("MS building for documents is not now available")
        return
    else:
        make_cmd = ['make', 'html']
    #  os.chdir(os.path.join(path, 'OneDQuantum/docs'))
    os.chdir(os.path.join(path, 'docs'))
    print("Building Documents")
    try:
        subprocess.check_call(make_cmd)
    except subprocess.CalledProcessError:
        print("Building documents failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Install ErwinJr2.')
    parser.add_argument('--msbuild',
                        help='Microsoft Visual Studio MSBuild directory. Use '
                        'Visual Studio compiler by providing this argument.',
                        type=str, nargs=1, default=None)
    parser.add_argument('--docs', help='Generate local documents.',
                        action='store_true', default=argparse.SUPPRESS)
    parser.add_argument('--nodocs', help='Do not generate local documents.',
                        dest='docs', action='store_false',
                        default=argparse.SUPPRESS)
    parser.add_argument('--shortcut', help='Generate shortcuts on desktop.',
                        action='store_true', default=argparse.SUPPRESS)
    parser.add_argument('--noshortcut',
                        help='Do not generate shortcuts on desktop.',
                        dest='shortcut', action='store_false',
                        default=argparse.SUPPRESS)

    args = parser.parse_args()

    currentPath = os.path.dirname(os.path.abspath(__file__))
    build_clib(currentPath, args.msbuild)
    if not hasattr(args, 'docs'):
        key = None
        while key not in ('y', 'n', ''):
            key = input("Build the documentation locally? Y/[N] ").lower()
        args.docs = (key == 'y')
    if args.docs:
        build_doc(currentPath)
    if not hasattr(args, 'shortcut'):
        key = None
        while key not in ('y', 'n', ''):
            key = input("Create shortcut on Desktop? [Y]/N ").lower()
        args.shortcut = not (key == 'n')
    if args.shortcut:
        try:
            create_shortcut(currentPath)
        except (NotImplementedError, CalledProcessError):
            print("Creat shortcut failed.")

# vim: ts=4 sw=4 sts=4 expandtab
