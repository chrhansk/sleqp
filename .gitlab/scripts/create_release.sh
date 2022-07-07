#!/usr/bin/env bash

# Fail the script when expansions fail
shopt -s failglob

set -x

VERSION_FILE="build/version"

VERSION=$(< ${VERSION_FILE})

TAG_NAME=v${VERSION}

if [ "$CI_COMMIT_BRANCH" != "releases/$VERSION" ]; then
  echo "Branch name $CI_COMMIT_BRANCH inconsistent with version $VERSION"
  exit 1
fi

# Check if release already exists

git tag | grep $VERSION

if [ "$?" == "0" ]; then
  echo "Duplicate release"
  exit 1
fi

set -e

PROJECT_URL="${CI_API_V4_URL}/projects/${CI_PROJECT_ID}"
TAG_URL="${PROJECT_URL}/repository/tags"
RELEASE_URL="${PROJECT_URL}/releases"
PACKAGE_URL="${PROJECT_URL}/packages/generic"

upload_file()
{
  file=$1
  shift
  group=$1

  filename=$(basename ${file})

  upload_url="${PACKAGE_URL}/${group}/${VERSION}/${filename}"

  curl --fail                                                  \
     --header "PRIVATE-TOKEN: ${SLEQP_RELEASE_TOKEN}"          \
     --upload-file ${file} ${upload_url} 1>&2

  link_name=$(cat <<EOF
{
  "name": "${filename}",
  "url": "${upload_url}",
  "filepath": "/${group}/${filename}",
  "link_type": "other"
}
EOF
)
  echo $link_name
}

# Create tags

curl --fail                                                    \
     --header "PRIVATE-TOKEN: ${SLEQP_RELEASE_TOKEN}"          \
     --request POST                                            \
     "${TAG_URL}?tag_name=${TAG_NAME}&ref=${CI_COMMIT_BRANCH}"

# Upload files to registry

link_names=""

for file in build/*.deb; do
  link_name=$(upload_file ${file} "debian")
  if [ "${link_names}" != "" ]; then
     link_names="${link_names}, ${link_name}"
  else
    link_names="${link_name}"
  fi
done

for file in bindings/python/dist/*.tar.*; do
  link_name=$(upload_file ${file} "python")
  if [ "${link_names}" != "" ]; then
     link_names="${link_names}, ${link_name}"
  else
    link_names="${link_name}"
  fi
done

for file in build/sleqp_octave*.tar.gz; do
  link_name=$(upload_file ${file} "octave_mex")
  if [ "${link_names}" != "" ]; then
     link_names="${link_names}, ${link_name}"
  else
    link_names="${link_name}"
  fi
done

body=$(cat <<EOF
{
  "name": "Release ${VERSION}",
  "tag_name": "${TAG_NAME}",
  "assets": {
    "links": [
      ${link_names}
    ]
  }
}
EOF
)

# Create release

curl --fail                                                    \
     --header "PRIVATE-TOKEN: ${SLEQP_RELEASE_TOKEN}"          \
     --header "Content-Type: application/json"                 \
     --data "${body}"                                          \
     "${RELEASE_URL}"
