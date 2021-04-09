"""
This file provides a detection of dark mode for MacOS
"""
# The code is a modified version of darkdetect
# (https://github.com/albertosottile/darkdetect),
# which is under the following licence

# Copyright (c) 2019, Alberto Sottile
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of "darkdetect" nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL "Alberto Sottile" BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# © 2020 GitHub, Inc.

import sys
import platform

isdark = False

if sys.platform == "darwin":
    from distutils.version import LooseVersion as V
    if V(platform.mac_ver()[0]) >= V("10.14"):
        import ctypes
        import ctypes.util

        appkit = ctypes.cdll.LoadLibrary(ctypes.util.find_library('AppKit'))
        objc = ctypes.cdll.LoadLibrary(ctypes.util.find_library('objc'))

        void_p = ctypes.c_void_p
        ull = ctypes.c_uint64

        objc.objc_getClass.restype = void_p
        objc.sel_registerName.restype = void_p
        objc.objc_msgSend.restype = void_p
        objc.objc_msgSend.argtypes = [void_p, void_p]

        msg = objc.objc_msgSend

        def _utf8(s):
            if not isinstance(s, bytes):
                s = s.encode('utf8')
            return s

        def n(name):
            return objc.sel_registerName(_utf8(name))

        def C(classname):
            return objc.objc_getClass(_utf8(classname))

        NSAutoreleasePool = objc.objc_getClass('NSAutoreleasePool')
        pool = msg(NSAutoreleasePool, n('alloc'))
        pool = msg(pool, n('init'))

        NSUserDefaults = C('NSUserDefaults')
        stdUserDef = msg(NSUserDefaults, n('standardUserDefaults'))

        NSString = C('NSString')

        key = msg(NSString, n("stringWithUTF8String:"), _utf8('AppleInterfaceStyle'))
        appearanceNS = msg(stdUserDef, n('stringForKey:'), void_p(key))
        appearanceC = msg(appearanceNS, n('UTF8String'))

        if appearanceC is not None:
            out = ctypes.string_at(appearanceC)
        else:
            out = None

        msg(pool, n('release'))

        if out is not None:
            isdark = out.decode('utf-8') == 'Dark'

    del V

del sys, platform
