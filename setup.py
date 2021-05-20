#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
from wheel.bdist_wheel import bdist_wheel, get_platform
import os
import subprocess
import warnings


_built = False


def build_clib():
    global _built
    if _built:
        return
    _built = True
    print("Building binary for ErwinJr2.")
    cwd = os.getcwd()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'ErwinJr2', 'OneDQuantum')
    os.chdir(path)
    # compile binary
    MSBuild = os.environ.get('MSBUILD')
    if MSBuild is not None:
        sln = 'OneDQuantum.sln'
        makeCMD = [MSBuild, sln, '/p:Configuration=Release']
        makeCMD_MP = [MSBuild, sln, '/p:Configuration=MP_Release']
    else:
        makeCMD = ['make']
        makeCMD_MP = ['make', 'MP']
    try:
        print('call ', makeCMD)
        subprocess.check_call(makeCMD)
    except subprocess.CalledProcessError:
        warnings.warn("Warning: Exit without compiling the C library. "
                      "Features are limited.")
        os.chdir(cwd)
        return
    try:
        print('call ', makeCMD_MP)
        subprocess.check_call(makeCMD_MP)
    except subprocess.CalledProcessError:
        warnings.warn("Warning: openMP not supported")
    os.chdir(cwd)


class EJBdistCMD(bdist_wheel):
    def run(self):
        # This is a hack s.t. the only way to let pip call build and use the
        # result
        build_clib()
        super().run()

    def finalize_options(self):
        super().finalize_options()
        self.plat_name_supplied = True
        self.plat_name = get_platform(self.bdist_dir)


class EJInstallCMD(install):
    def run(self):
        build_clib()
        super().run()


class EJDevelopCMD(develop):
    def run(self):
        build_clib()
        super().run()


long_description = """
This is a Quantum Cascade Laser (QCL) modeling and design software produced
at Princeton University, Gmachl group.

See https://erwinjr2.readthedocs.io/ for details.
"""

setup(
    name='ErwinJr2',
    # Important! keep Version consistent for major and minor with
    # versionAndName.py as well as the release tag on GitHub
    version='2.2.3',
    author='Ming Lyu',
    author_email='minglyu@princeton.edu',
    license="GPL-3.0",
    url='https://princetonuniversity.github.io/ErwinJr2',
    description='A Quantum Cascade Laser Design Tool',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    package_data={
        'ErwinJr2': [
            'images/*.png', 'images/*.ico', 'images/*.icns',
            'example/PQLiu.json',
            'Info.plist'
        ],
        'ErwinJr2.OneDQuantum': [
            'Makefile', '*.sln', '*.vcxproj',
            '*.c', '*.h', 'fftautocorr/*.c', 'fftautocorr/*.h',
            '*.so', '*.dll', '*.dylib'
        ]
    },
    cmdclass={
        'bdist_wheel': EJBdistCMD,
        'install': EJInstallCMD,
        'develop': EJDevelopCMD
    },
    python_requires='>=3.6',
    install_requires=['numpy>=1.12', 'scipy>=0.18'],
    extras_require={
        'GUI': [
            'pyqt5',
            'matplotlib>=3.3',
            "winshell;platform_system=='Windows'"
        ],
        'dev': [
            # Required by documentation
            'sphinx>=2',
            'sphinx_rtd_theme',
            'sphinxcontrib-bibtex>=2.0.0',
            'breathe'
        ]
    },
    entry_points={
        'console_scripts': [
            'ErwinJr = ErwinJr2.__main__:__main__ [GUI]',
            'ErwinJr-genshortcut = ErwinJr2.genshortcut:create_shortcut [GUI]'
        ]
    }
)
