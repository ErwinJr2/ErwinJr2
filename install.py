#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# TODO: add dependence check
import argparse

import os
import sys
import subprocess
from subprocess import CalledProcessError
if sys.platform.startswith('win'):
    try:
        import winshell  # ignore: unresolved-import
    except ModuleNotFoundError:
        print("winshell not installed. Cannot create shortcut!")

verbose = False


def build_clib(path, MSBuild=None):
    if MSBuild is None:
        make_cmd = ['make']
        makemp_cmd = ['make', 'MP']
    else:
        make_cmd = [MSBuild, 'OneDQuantum.sln', '/p:Configuration=Release']
        makemp_cmd = [MSBuild, '1DSchrodinger.vcxproj',
                      '/p:Configuration=MP_Release']
        # print('running cummand: ')
        # print(' '.join(make_cmd))
    os.chdir(os.path.join(path, 'OneDQuantum'))
    print("Building C Lib")
    result = subprocess.run(make_cmd, check=True)
    print(result.stdout)
    print(result.stderr)
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


def create_shortcut(path):
    os.chdir(path)
    if sys.platform.startswith('darwin'):
        shortcut_path = os.path.join(os.environ["HOME"],
                                     "Desktop", "ErwinJr2.app", "Contents")
        subprocess.check_call(['mkdir', '-p', shortcut_path])
        os.chdir(shortcut_path)
        subprocess.check_call(['mkdir', 'MacOS'])
        subprocess.check_call(['mkdir', 'Resources'])
        subprocess.check_call(['cp', os.path.join(path, 'Info.plist'), './'])
        subprocess.check_call([
            'cp', os.path.join(path, 'images', 'EJicns.icns'), 'Resources/'])
        shellcode = """#!/bin/bash\n"""
        shellcode += "cd %s\n" % path
        shellcode += "%s ErwinJr.pyw\n" % sys.executable
        shellfile = 'MacOS/ErwinJr2'
        with open(shellfile, 'w') as f:
            f.write(shellcode)
        subprocess.call(['chmod', '+x', shellfile])
    if sys.platform.startswith('win'):
        shortcut_path = os.path.join(os.path.expanduser('~'),
                                     'Desktop', 'ErwinJr2.lnk')
        batcode = """@SET "PATH=%s\"\n""" % os.environ['PATH']
        batcode += "@echo off\n"
        batcode += "cd %~dp0\n"
        batcode += "start pythonw ErwinJr.pyw %1\n"
        batfile = "ErwinJr2.bat"
        with open(batfile, 'w') as f:
            f.write(batcode)
        with winshell.shortcut(shortcut_path) as link:
            link.path = os.path.join(path, batfile)
            link.description = "Shortcut to ErwinJr2"
            link.icon_location = (os.path.join(path, 'images', 'EJico.ico'), 0)
            link.working_directory = path
            link.show_cmd = 'min'
    else:
        raise NotImplementedError("Operating system not supported.")


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
    parser.add_argument('--verbose', '-v', action='count', default=0)

    args = parser.parse_args()
    verbose = args.verbose > 0

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
