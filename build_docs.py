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
        "name": "latest",
        "version": "latest",
        "url": "https://sandialabs.github.io/sansmic/",
        "preferred": True,
    }
]

for tag in tags:
    versions.append(
        dict(
            name=tag,
            version=tag,
            url="https://sandialabs.github.io/sansmic/" + tag + "/",
            preferred=False,
        )
    )

with open(
    os.path.abspath(os.path.join(".", "docs", "_static", "switcher.json")), "w"
) as fswitch:
    json.dump(versions, fswitch)

my_env["SANSMIC_SPHINX_VERSION"] = "latest"

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

for tag in tags:
    os.environ["VERSION_INFO"] = repr(tag)
    my_env["SANSMIC_SPHINX_VERSION"] = tag
    subprocess.run("git checkout " + tag, shell=True)
    subprocess.run(
        "git checkout main -- docs/conf.py",
        shell=True,
    )
    subprocess.run(
        "pip install --no-deps --force-reinstall --no-input --no-cache-dir -e .",
        shell=True,
    )
    files = glob.glob("./docs/apidocs/*.rst")
    for f in files:
        os.remove(f)
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
