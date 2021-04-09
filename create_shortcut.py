#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
if sys.platform.startswith('win'):
    try:
        import winshell  # type: ignore # ignore: unresolved-import
    except ModuleNotFoundError:
        print("winshell not installed. Cannot create shortcut!")


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
    currentPath = os.path.dirname(os.path.abspath(__file__))
    create_shortcut(currentPath)

# vim: ts=4 sw=4 sts=4 expandtab
