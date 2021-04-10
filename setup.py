#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import sys
import subprocess


class EJBuildBinaryCMD(install):
    def run(self):
        install.run(self)
        cwd = os.getcwd()
        path = os.path.dirname(os.path.abspath(__file__))
        print("Building binary for ErwinJr2.")
        # compile binary
        MSBuild = os.environ.get('MSBUILD')
        if MSBuild is not None:
            make_cmd = [MSBuild, 'OneDQuantum.sln', '/p:Configuration=Release']
            makemp_cmd = [MSBuild, '1DSchrodinger.vcxproj',
                          '/p:Configuration=MP_Release']
        else:
            make_cmd = ['make']
            makemp_cmd = ['make', 'MP']
        os.chdir(os.path.join(path, 'ErwinJr2/OneDQuantum'))
        print("Building C Lib")
        try:
            subprocess.check_call(make_cmd)
        except subprocess.CalledProcessError:
            print("Exit without compiling the C library. "
                  "Features are limited.")
            os.chdir(cwd)
            return
        try:
            subprocess.check_call(makemp_cmd)
        except subprocess.CalledProcessError:
            print("openMP not supported")
        os.chdir(cwd)


setup(
    name='ErwinJr2',
    version='2.0.0',
    author='Ming Lyu',
    author_email='minglyu@princeton.edu',
    license="GPL-3.0",
    url='https://princetonuniversity.github.io/ErwinJr2',
    description='A Quantum Cascade Laser Design Tool',
    packages=find_packages(),
    package_data={
        'ErwinJr2': [
            'images/*.png', 'images/*.ico', 'images/*.icns', 'example/*',
            'Info.plist'
        ],
        'ErwinJr2.OneDQuantum': ['*.so', '*.dll', '*.dylib']
    },
    cmdclass={
        'install': EJBuildBinaryCMD,
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
