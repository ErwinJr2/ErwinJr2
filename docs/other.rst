Version Control and Automatic Testing
======================================
We use ``git`` as the version control tool, and the git repository is hosted on 
`github <https://github.com/PrincetonUniversity/OneDQ/>`_.

``master`` branch is meant to be a stable branch. 
The development happens in ``dev`` branch.
``doc`` branch is for design document and is no longer maintained. New documents 
are embedded within the code using `sphinx` for Python and `doxygen` + `breathe` 
for C. 

We also have automatic testing system using `travis-ci <http://travis-ci.com/>`_ 
service. The testing code is under ``test/`` directory and ``.travis.yml`` script. 
