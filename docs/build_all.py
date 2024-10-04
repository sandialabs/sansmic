# coding: utf-8

"""Get all tagged releases and build the documentation.
"""

import subprocess
import json
import os.path

git_tag = subprocess.run(
    ["git", "tag", "--list", '"v*.*.*"', '--sort="version:refname"'],
    capture_output=True,
    text=True,
)
tags = git_tag.stdout.splitlines()
versions = [
    {
        "name": "latest",
        "version": "latest",
        "url": "https://sandialabs.github.io/sansmic/",
        "preferred": True,
    }
]

subprocess.run(["sphinx-build", "-b", "html", ".", os.path.join("_build", "html")])

for tag in tags:
    versions.append(
        dict(
            name=tag,
            version=tag,
            url="https://sandialabs.github.io/sansmic/" + tag + "/",
            preferred=False,
        )
    )
    os.environ["VERSION_INFO"] = repr(tag)
    subprocess.run(["git", "checkout", tag])
    subprocess.run(
        ["sphinx-build", "-b", "html", ".", os.path.join("_build", "html", tag)]
    )

with open(os.path.join("_build", "html", "_static", "switcher.json"), "w") as fswitch:
    json.dump(versions, fswitch)
