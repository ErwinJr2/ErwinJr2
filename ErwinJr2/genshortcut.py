#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
if sys.platform.startswith('win'):
    import winshell  # type: ignore # ignore: unresolved-import


def macos_shortcut(basePath, shortcut_path=None):
    os.chdir(basePath)
    if shortcut_path is None:
        shortcut_path = os.path.join(os.environ["HOME"], "Desktop")
    contents = os.path.join(shortcut_path, "ErwinJr2.app", "Contents")
    if os.path.exists(contents):
        subprocess.call(['rm', '-r', contents])
    subprocess.check_call(['mkdir', '-p', contents])
    os.chdir(contents)
    subprocess.check_call(['mkdir', 'MacOS'])
    subprocess.check_call(['mkdir', 'Resources'])
    subprocess.check_call(['cp', os.path.join(basePath, 'Info.plist'), './'])
    subprocess.check_call([
        'cp', os.path.join(basePath, 'images', 'EJicns.icns'), 'Resources/'])
    shellcode = "#!/bin/bash\n"
    if os.environ.get('PYTHONPATH'):
        shellcode += "export PYTHONPATH=%s\n" % os.environ['PYTHONPATH']
    shellcode += "%s -m ErwinJr2\n" % sys.executable
    shellfile = 'MacOS/ErwinJr2'
    with open(shellfile, 'w') as f:
        f.write(shellcode)
    subprocess.call(['chmod', '+x', shellfile])


def win_shortcut(basePath, shortcut_path):
    os.chdir(basePath)
    if shortcut_path is None:
        shortcut_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    link_path = os.path.join(shortcut_path, 'ErwinJr2.lnk')
    batcode = """@SET "PATH=%s\"\n""" % os.environ['PATH']
    if os.environ.get('PYTHONPATH'):
        batcode += """@SET "PYTHONPATH=%s\"\n""" % os.environ['PYTHONPATH']
    # batcode += "@echo off\n"
    # batcode += "cd %~dp0\n"
    batcode += "start pythonw -m ErwinJr2 %1\n"
    batfile = "ErwinJr2.bat"
    with open(batfile, 'w') as f:
        f.write(batcode)
    with winshell.shortcut(link_path) as link:
        link.path = os.path.join(basePath, batfile)
        link.description = "Shortcut to ErwinJr2"
        link.icon_location = (os.path.join(basePath, 'images', 'EJico.ico'), 0)
        link.working_directory = basePath
        link.show_cmd = 'min'


linux_app = """#!/usr/bin/env xdg-open
[Desktop Entry]
Type=Application
Name=ErwinJr2
Comment=ErwinJr2 Launcher
Exec=python3 -m ErwinJr2
Icon={icon}
Terminal=false
"""


def linux_shortcut(basePath, shortcut_path):
    os.chdir(basePath)
    if shortcut_path is None:
        shortcut_path = os.path.join(
            os.path.expanduser('~'), '.local', 'share', 'applications')
    app_path = os.path.join(shortcut_path, 'ErwinJr2.desktop')
    icon_path = os.path.join(basePath, 'images', 'EJpng256.png')
    app = linux_app.format(icon=icon_path)
    with open(app_path, 'w') as f:
        f.write(app)
    subprocess.call(['chmod', '+x', app_path])
    try:
        desktop_path = subprocess.check_output(['xdg-user-dir', 'DESKTOP'])
        desktop_path = desktop_path.decode().split()[0]
        assert(os.path.isdir(desktop_path))
        subprocess.call(['cp', app_path, desktop_path])
    except (FileNotFoundError, IndexError, AssertionError):
        print('Desktop directory not found.')


def create_shortcut(shortcut_path=None):
    basePath = os.path.dirname(os.path.abspath(__file__))
    if sys.platform.startswith('darwin'):
        macos_shortcut(basePath, shortcut_path)
    elif sys.platform.startswith('win'):
        win_shortcut(basePath, shortcut_path)
    elif sys.platform.startswith('linux'):
        linux_shortcut(basePath, shortcut_path)
    else:
        raise NotImplementedError("Operating system not supported.")


if __name__ == "__main__":
    create_shortcut()
