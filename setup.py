#!/usr/bin/env python
## -*- encoding: utf-8 -*-

import os
import re
import sys
from setuptools import setup
from codecs import open  # To open the README file with proper encoding
from setuptools.command.test import test as TestCommand  # for tests
from setuptools.extension import Extension
from sage.env import sage_include_directories
from Cython.Build import cythonize

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("sage -t --force-lib pyhdme")
        if errno != 0:
            sys.exit(1)


cythonize_dir = "build"

path = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(path, "pyhdme/lib")
data_path = os.path.join(lib_path, "hdme_data")
allfiles_in_lib = [
    os.path.relpath(os.path.join(dp, f), path)
    for dp, dn, fn in os.walk(os.path.expanduser(lib_path))
    for f in fn
]

hdme_sources = [
    elt
    for elt in allfiles_in_lib
    if elt.endswith(".c") and
    "test/" not in elt and
    "time/" not in elt and
    "examples/" not in elt and
    "programs/" not in elt
]

hdme_data_files = [
    os.path.relpath(elt, data_path)
    for elt in allfiles_in_lib
    if elt.startswith(data_path) and not elt.endswith(".c")
]






def patch_hdme_data_read(data_files):
    def to_C(filename):
        value = readfile(os.path.join(data_path, filename))
        hex_value = [format(c, '#04x') for c in value.encode('ascii')] + ['0x00'] # null terminated
        varname = filename.replace('/', '_')
        return f'/*\n{value}\n*/\nstatic const char {varname}[] = {{{", ".join(hex_value)}}};\n'

    headers = """
#include <stdint.h>
#include <string.h>

#include "hdme_data.h"
"""
    C_lookup = "{\n" + ",\n".join([
        f'{{"{elt}", {elt.replace("/", "_")}}}'
        for elt in data_files]) + "\n}"
    lookup_table = """
typedef struct { char *key; const char *val; } pair;
static pair lookuptable[] = %s;
#define NKEYS (sizeof(lookuptable)/sizeof(pair))

char* get_string_from_key(char *dest, const char *key)
{
    size_t i; /* GCC does not allow variable declarations in for loop initializers before C99 */
    for(i=0; i < NKEYS; ++i) {
        pair *item = lookuptable + i;
        if (strcmp(item->key, key) == 0)
            return strcpy(dest, item->val);
    }
    return NULL;
}

""" % (C_lookup,)

    filename = "pyhdme/lib/hdme_data/hdme_data_read.c"
    newfilename = "pyhdme/lib/hdme_data/hdme_data_read_static.c"
    oldfilename = filename + ".old"
    if not os.path.exists(oldfilename):
        os.rename(filename, oldfilename)
    source = readfile(oldfilename)
    # comment out lines starting with
    for elt in ['#include', 'char filename', 'FILE', 'file', 'flint_sprintf', 'fclose']:
        source = re.sub(fr'^(\s+)({elt}.*)$', r'\1/*\2*/', source, count=1, flags=re.MULTILINE)
    source = re.sub(r'^(\s+)(success.*)$', r'\1/*\2*/\n\1success = get_string_from_key(str, name);\n',
                    source,
                    count=1,
                    flags=re.MULTILINE)
    source = re.sub(fr'^(.*)(filename)(\);)$', r'\1name\3', source, flags=re.MULTILINE)
    source = "\n".join([headers, "".join(map(to_C, data_files)), lookup_table, source])

    with open(newfilename, 'w') as W:
        W.write(source)
    return filename, newfilename

def undo_patch_hdme_data_read():
    filename = "pyhdme/lib/hdme_data/hdme_data_read.c"
    newfilename = "pyhdme/lib/hdme_data/hdme_data_read_static.c"
    oldfilename = filename + ".old"
    if os.path.exists(oldfilename):
        os.rename(oldfilename, filename)
    if os.path.exists(newfilename):
        os.remove(newfilename)


old, new = patch_hdme_data_read(hdme_data_files)
if old in hdme_sources:
    hdme_sources.remove(old)
if new not in hdme_sources:
    hdme_sources.append(new)


pyhdme = Extension(
    "pyhdme.hdme",
    language="c",
    sources=[
        "pyhdme/hdme.pyx",
    ]
    + hdme_sources,
    libraries=["arb", "flint", "mpfr", "gmp", "pthread", "m"],
    include_dirs=sage_include_directories() + ["pyhdme/lib/"],
)

setup(
    name="pyhdme",
    author="Edgar Costa",
    author_email="edgarcosta@math.dartmouth.edu",
    url="https://github.com/edgarcosta/pyhdme",
    license="GNU General Public License, version 3",
    description="Wrapper for C library for evaluating higher-dimensional modular equations",
    long_description=readfile("README.md"),  # get the long description from the README
    version=readfile("VERSION"),  # the VERSION file is shared with the documentation
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: GNU General Public License v2 or v3",
        "Programming Language :: Python :: 3.7",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords="sagemath hdme",
    setup_requires=[
        "cython",
        "sagemath",
    ],  # currently useless, see https://www.python.org/dev/peps/pep-0518/
    install_requires=["cython", "sagemath", "sphinx"],
    packages=["pyhdme"],
    include_package_data=False,
    ext_modules=cythonize([pyhdme], language="c"),
    cmdclass={"test": SageTest}  # adding a special setup command for tests
    # ext_modules = extensions,
    # cmdclass = {'test': SageTest, 'build_ext': Cython.Build.build_ext} # adding a special setup command for tests and build_ext
)

undo_patch_hdme_data_read()
