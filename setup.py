#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# import wheel.bdist_wheel as bdist_wheel
# from wheel.bdist_wheel import bdist_wheel, get_platform
import os
import subprocess
import warnings

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.dist import Distribution


def build_clib():
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


class EJDevelopCMD(develop):
    def run(self):
        build_clib()
        super().run()


class EJBuildCMD(build_py):
    def run(self):
        build_clib()
        super().run()


class EJDistribution(Distribution):
    def has_ext_modules(self):
        super().has_ext_modules()
        return True


setup(
    cmdclass={
        'develop': EJDevelopCMD,
        'build_py': EJBuildCMD,
    },
    distclass=EJDistribution,
)
