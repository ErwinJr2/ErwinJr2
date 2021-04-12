echo Depolying $TRAVIS_BRANCH $TRAVIS_TAG on $TRAVIS_OS_NAME
if [ $TRAVIS_BRANCH = 'dev' ]; then
    $PYEXE setup.py $1
    $PYEXE -m twine upload --repository testpypi dist/* \
        --username __token__ --password "$PYPITEST_PASSWORD" \
        --skip-existing
elif [ ! -z $TRAVIS_TAG ]; then
    $PYEXE setup.py $1
    $PYEXE -m twine upload dist/* \
        --username __token__ --password "$PYPI_PASSWORD" \
        --skip-existing
else
    echo skip deploy
fi
