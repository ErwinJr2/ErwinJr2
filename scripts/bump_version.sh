#!/bin/bash

VERSION_FILE=ErwinJr2/version.py

set -e -x
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$CURRENT_BRANCH" != "dev" ]; then
    echo "You must be on the dev branch to run this script"
    exit 1
fi

DEV_VERSION=$(grep -oP 'VERSION = "\K[^"]+' $VERSION_FILE)
VERSION_MAJOR="${DEV_VERSION%%.*}"
VERSION_MINOR_PATCH="${DEV_VERSION#*.}"
VERSION_MINOR="${VERSION_MINOR_PATCH%%.*}"
VERSION_PATCH_PRE_RELEASE="${VERSION_MINOR_PATCH#*.}"
VERSION_PATCH="${VERSION_PATCH_PRE_RELEASE%%-dev}"
VERSION="$VERSION_MAJOR.$VERSION_MINOR.$VERSION_PATCH"
NEW_DEV_VERSION=$VERSION_MAJOR.$VERSION_MINOR.$((VERSION_PATCH+1))-dev

sed -i "s/VERSION = \"$DEV_VERSION\"/VERSION = \"$VERSION\"/" $VERSION_FILE
git add $VERSION_FILE
git commit -m "Bump version to $VERSION"
git checkout master
git merge dev
git tag -a v$VERSION -m "Release v$VERSION"
git push --atomic origin master v$VERSION

git branch -D dev
git checkout -b dev
sed -i "s/VERSION = \"$VERSION\"/VERSION = \"$NEW_DEV_VERSION\"/" $VERSION_FILE
git commit -am "Bump version to $NEW_DEV_VERSION"
git push --set-upstream -f origin dev
