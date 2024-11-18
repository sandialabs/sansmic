# coding: utf-8

"""Get all tagged releases and build the documentation.
"""

import subprocess
import json
import os
import glob

my_env = os.environ.copy()

git_tag = subprocess.run(
    " ".join(["git", "tag", "--list", "v*.*.*", "--sort=-version:refname"]),
    capture_output=True,
    text=True,
    shell=True,
)

tags = git_tag.stdout.splitlines()
versions = [
    {
        "name": "main",
        "version": "main",
        "url": "https://sandialabs.github.io/sansmic/",
        "preferred": False,
    }
]

found_stable = False
latest = False
for tag in tags:
    if not found_stable and not "-" in tag:
        latest = True
        found_stable = tag
    else:
        latest = False
    versions.append(
        dict(
            name=tag,
            version=tag,
            url="https://sandialabs.github.io/sansmic/" + tag + "/",
            preferred=latest,
        )
    )

with open(
    os.path.abspath(os.path.join(".", "docs", "_static", "switcher.json")),
    "w",
) as fswitch:
    json.dump(versions, fswitch)
