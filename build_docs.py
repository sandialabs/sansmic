# coding: utf-8

"""Get all tagged releases and build the documentation.
"""

import subprocess
import json
import os
import glob

git_tag = subprocess.run(
    " ".join(["git", "tag", "--list", "v*.*.*", "--sort=version:refname"]),
    capture_output=True,
    text=True,
    shell=True,
)

my_env = os.environ.copy()

tags = git_tag.stdout.splitlines()
versions = [
    {
        "name": "main",
        "version": "main",
        "url": "https://sandialabs.github.io/sansmic/main/",
        "preferred": False,
    }
]

found_stable = False
latest = False
for tag in tags:
    if not found_stable and not "-" in tag:
        latest = True
        found_stable = tag
        versions.append(
            dict(
                name=latest,
                version=tag,
                url="https://sandialabs.github.io/sansmic/",
                preferred=True,
            )
        )
    else:
        latest = False
    versions.append(
        dict(
            name=tag,
            version=tag,
            url="https://sandialabs.github.io/sansmic/" + tag + "/",
            preferred=False,
        )
    )
if not found_stable:
    found_stable = "main"

with open(
    os.path.abspath(os.path.join(".", "docs", "_static", "switcher.json")), "w"
) as fswitch:
    json.dump(versions, fswitch)

tag = found_stable
os.environ["VERSION_INFO"] = repr(tag)
my_env["SANSMIC_SPHINX_VERSION"] = tag
files = glob.glob("./docs/apidocs/*.rst")
for f in files:
    os.remove(f)
subprocess.run("git checkout " + tag, shell=True)
subprocess.run(
    " ".join(
        [
            "sphinx-build",
            "-b",
            "html",
            "docs/",
            "docs/_build/html",
        ]
    ),
    shell=True,
    env=my_env,
)

os.remove(os.path.abspath(os.path.join(".", "docs", "_static", "switcher.json")))

tags.append("main")

for tag in tags:
    os.environ["VERSION_INFO"] = repr(tag)
    my_env["SANSMIC_SPHINX_VERSION"] = tag
    files = glob.glob("./docs/apidocs/*.rst")
    for f in files:
        os.remove(f)
    subprocess.run("git checkout " + tag, shell=True)
    subprocess.run(
        " ".join(
            [
                "sphinx-build",
                "-b",
                "html",
                "docs/",
                "docs/_build/html/" + tag,
            ]
        ),
        shell=True,
        env=my_env,
    )
