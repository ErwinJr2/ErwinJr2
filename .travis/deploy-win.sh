python setup.py bdist_wheel
python -m twine upload --repository testpypi dist/* --username __token__ --password "$PYPI_PASSWORD"