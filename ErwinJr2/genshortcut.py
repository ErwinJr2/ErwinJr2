#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""Scripts to generate shortcuts for ErwinJr2."""

import os
import subprocess
import sys

if sys.platform.startswith("win"):
    import winshell  # type: ignore # ignore: unresolved-import


def macos_shortcut(proj_root, shortcut_path=None):
    os.chdir(proj_root)
    if shortcut_path is None:
        shortcut_path = os.path.join(os.environ["HOME"], "Desktop")
    contents = os.path.join(shortcut_path, "ErwinJr2.app", "Contents")
    if os.path.exists(contents):
        subprocess.call(["rm", "-r", contents])
    subprocess.check_call(["mkdir", "-p", contents])
    os.chdir(contents)
    subprocess.check_call(["mkdir", "MacOS"])
    subprocess.check_call(["mkdir", "Resources"])
    subprocess.check_call(["cp", os.path.join(proj_root, "Info.plist"), "./"])
    subprocess.check_call(
        ["cp", os.path.join(proj_root, "images", "EJicns.icns"), "Resources/"]
    )
    shellcode = "#!/bin/bash\n"
    if python_path := os.environ.get("PYTHONPATH"):
        shellcode += f"export PYTHONPATH={python_path}\n"
    shellcode += f"{sys.executable} -m ErwinJr2\n"
    shellfile = "MacOS/ErwinJr2"
    with open(shellfile, "w") as f:
        f.write(shellcode)
    subprocess.call(["chmod", "+x", shellfile])


def win_shortcut(proj_root, shortcut_path):
    assert sys.platform.startswith("win"), "This function is for Windows only."
    os.chdir(proj_root)
    if shortcut_path is None:
        shortcut_path = os.path.join(os.path.expanduser("~"), "Desktop")
    link_path = os.path.join(shortcut_path, "ErwinJr2.lnk")
    batcode = f"""@SET "PATH={os.environ["PATH"]}\"\n"""
    if python_path := os.environ.get("PYTHONPATH"):
        batcode += f"""@SET "PYTHONPATH={python_path}\"\n"""
    batcode += "start pythonw -m ErwinJr2 %1\n"
    batfile = "ErwinJr2.bat"
    with open(batfile, "w") as f:
        f.write(batcode)
    with winshell.shortcut(link_path) as link:  # pylint: disable=E0606
        link.path = os.path.join(proj_root, batfile)
        link.description = "Shortcut to ErwinJr2"
        link.icon_location = (os.path.join(proj_root, "images", "EJico.ico"), 0)
        link.working_directory = proj_root
        link.show_cmd = "min"


linux_app = """#!/usr/bin/env xdg-open
[Desktop Entry]
Type=Application
Name=ErwinJr2
Comment=ErwinJr2 Launcher
Exec=python3 -m ErwinJr2
Icon={icon}
Terminal=false
"""


def linux_shortcut(proj_root, shortcut_path):
    os.chdir(proj_root)
    if shortcut_path is None:
        shortcut_path = os.path.join(
            os.path.expanduser("~"), ".local", "share", "applications"
        )
    app_path = os.path.join(shortcut_path, "ErwinJr2.desktop")
    icon_path = os.path.join(proj_root, "images", "EJpng256.png")
    app = linux_app.format(icon=icon_path)
    with open(app_path, "w") as f:
        f.write(app)
    subprocess.call(["chmod", "+x", app_path])
    try:
        desktop_path = subprocess.check_output(["xdg-user-dir", "DESKTOP"])
        desktop_path = desktop_path.decode().split()[0]
        assert os.path.isdir(desktop_path)
        subprocess.call(["cp", app_path, desktop_path])
    except (FileNotFoundError, IndexError, AssertionError):
        print("Desktop directory not found.")


def create_shortcut(shortcut_path=None):
    proj_root = os.path.dirname(os.path.abspath(__file__))
    if sys.platform.startswith("darwin"):
        macos_shortcut(proj_root, shortcut_path)
    elif sys.platform.startswith("win"):
        win_shortcut(proj_root, shortcut_path)
    elif sys.platform.startswith("linux"):
        linux_shortcut(proj_root, shortcut_path)
    else:
        raise NotImplementedError(f"OS not supported: {sys.platform}")


if __name__ == "__main__":
    create_shortcut()
