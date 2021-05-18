Version Control, Continuous Integration and Deployment
=======================================================

For the use of ``git`` version control, there are two major branches,
``master`` and ``dev``. The ``master`` branch is the default
branch when a user downloads the source code from GitHub. It is recommended to
keep this branch a relatively stable branch and only merge in new codes
when they are reliable. The ``dev`` branch is intended to be the
developing branch for a single maintainer.

For the version number, we suggest the three-part version number convention:
[major version.minor version.patch]. For example, version 2.1.0 means major
version 2, minor version 1, and patch 0. Major version should increase when
there are large updates or a change of old APIs; minor version should increase
when there are new features; patch number increases when it's a bug fix.
To tag a version, it is done on the GitHub repository webpage, so that new
versions will trigger a release and the following deployment.

The software requires optionally a compiled C library. For users installing the
software directly from `The Python Package Index <https://pypi.org/project/ErwinJr2/>`_ (PyPI) via ``pip``
command, Python will first check if a compiled binary file (called ``wheel``)
for the C library exists on PyPI and download the binary if so.
Otherwise, it will try to compile the C library locally from the source code.
We call uploading the source code and the compiled binary for Windows and
macOS ``deployment``.

As defined in ``.travis.yml``, for each commit pushed to GitHub, Travis-CI
will first run the tests under different platforms, notify the repository owner
if tests fail, and try to deploy the code to PyPI.
The branches also influence the behavior of the deployment. The current
setting is that, for each new commit to ``dev`` branch, Travis-CI will
try to deploy the software on `test.pypi <test.pypi.org/project/ErwinJr2/>`_ if there is
not a duplicate version already, which usually means the version number defined
in ``setup.py`` is updated.
For the ``master``, Travis-CI will not try to deploy the software
to the official PyPI unless it's tagged on GitHub
as a new release.

Apart from Travis-CI, ``readthedocs`` also monitors new commits to the
GitHub. It will create the online document you are now reading for any commits to the ``master``
and the ``dev`` branch, creating ``latest`` and ``dev`` version of the
document. For any tagged version it will create a ``stable`` version of the
document.
