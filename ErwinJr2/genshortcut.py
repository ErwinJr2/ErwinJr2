#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
if sys.platform.startswith('win'):
    import winshell  # type: ignore # ignore: unresolved-import


def create_shortcut(shortcut_path=None):
    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    if sys.platform.startswith('darwin'):
        if shortcut_path is None:
            shortcut_path = os.path.join(os.environ["HOME"], "Desktop")
        contents = os.path.join(shortcut_path, "ErwinJr2.app", "Contents")
        if os.path.exists(contents):
            subprocess.call(['rm', '-r', contents])
        subprocess.check_call(['mkdir', '-p', contents])
        os.chdir(contents)
        subprocess.check_call(['mkdir', 'MacOS'])
        subprocess.check_call(['mkdir', 'Resources'])
        subprocess.check_call(['cp', os.path.join(path, 'Info.plist'), './'])
        subprocess.check_call([
            'cp', os.path.join(path, 'images', 'EJicns.icns'), 'Resources/'])
        shellcode = """#!/bin/bash\n"""
        shellcode += "%s -m ErwinJr2\n" % sys.executable
        shellfile = 'MacOS/ErwinJr2'
        with open(shellfile, 'w') as f:
            f.write(shellcode)
        subprocess.call(['chmod', '+x', shellfile])
    elif sys.platform.startswith('win'):
        if shortcut_path is None:
            shortcut_path = os.path.join(os.path.expanduser('~'), 'Desktop')
        link_path = os.path.join(shortcut_path, 'ErwinJr2.lnk')
        batcode = """@SET "PATH=%s\"\n""" % os.environ['PATH']
        batcode += "@echo off\n"
        batcode += "cd %~dp0\n"
        batcode += "start pythonw ErwinJr.py %1\n"
        batfile = "ErwinJr2.bat"
        with open(batfile, 'w') as f:
            f.write(batcode)
        with winshell.shortcut(link_path) as link:
            link.path = os.path.join(path, batfile)
            link.description = "Shortcut to ErwinJr2"
            link.icon_location = (os.path.join(path, 'images', 'EJico.ico'), 0)
            link.working_directory = path
            link.show_cmd = 'min'
    else:
        raise NotImplementedError("Operating system not supported.")


if __name__ == "__main__":
    create_shortcut()

# vim: ts=4 sw=4 sts=4 expandtab
