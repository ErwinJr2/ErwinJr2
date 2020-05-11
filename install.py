#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# TODO: add dependence check

import os
import sys
import subprocess
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


def create_shortcut(path):
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
    else:
        raise NotImplementedError("Operating system not supported.")


if __name__ == "__main__":
    MSBuild = None
    path = os.path.dirname(os.path.abspath(__file__))
    docs = None
    shortcut = None
    for opt in sys.argv[1:]:
        if opt.lower().startswith("--msbuild="):
            MSBuild = opt[10:]
        elif opt.lower().startswith("--nodoc"):
            docs = False
        elif opt.lower().startswith("--doc"):
            docs = True
        elif opt.lower().startswith("--shortcut"):
            shortcut = True
        elif opt.lower().startswith("--noshortcut"):
            shortcut = False
        else:
            print("Unknown option %s" % opt)
    build_clib(path, MSBuild)
    if docs is None:
        key = None
        while key not in ('y', 'n', ''):
            key = input("Build the documentation locally? Y/[N] ").lower()
        docs = (key == 'y')
    if docs:
        build_doc(path, MSBuild)
    if shortcut is None:
        key = None
        while key not in ('y', 'n', ''):
            key = input("Create shortcut on Desktop? [Y]/N ").lower()
        shortcut = not (key == 'n')
    if shortcut:
        try:
            create_shortcut(path)
        except (NotImplementedError, CalledProcessError):
            print("Creat shortcut failed.")

# vim: ts=4 sw=4 sts=4 expandtab
