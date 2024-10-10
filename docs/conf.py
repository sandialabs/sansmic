# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

import os
import sys
import glob

sys.path.insert(0, os.path.abspath(os.path.join("..", "src", "python")))
import sansmic

doxygen_installed = (
    os.system("doxygen --version") == 0
)  # doxygen must be installed for C++ docs


##############################################################################
#    Sphinx core options                                                     #
##############################################################################

# -- Project information -----------------------------------------------------
project = "sansmic"
copyright = "2024 National Technology and Engineering Solutions of Sandia, LLC. Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software."
author = "See AUTHORS.md"
version = os.environ.get("SANSMIC_SPHINX_VERSION", sansmic.__version__)
release = version

# -- Extensions to load -------------------------------------------------------
extensions = [
    "breathe",
    "exhale",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "sphinx_design",
    "sphinx_click",
    "sphinxcontrib.bibtex",
]


# -- Callout numbering -------------------------------------------------------
numfig = True
numfig_format = {"figure": "Figure %s", "table": "Table %s", "code-block": "Listing %s"}


# -- Internationalization ----------------------------------------------------
language = "en"


# -- Markup options ----------------------------------------------------------


# -- Options for Maths -------------------------------------------------------


# -- Options for object signatures
add_function_parentheses = True
toc_object_entries = True
toc_object_entries_show_parents = "hide"


# -- Options for source files ------------------------------------------------

# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "doxygen"]

# The master toctree document.
master_doc = "index"

# source_suffix = ['.rst', '.md']
source_suffix = {".rst": "restructuredtext"}


# -- Options for templating --------------------------------------------------

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]


##############################################################################
#    Sphinx domain options                                                   #
##############################################################################

# -- Options for the Python domain -------------------------------------------
add_module_names = False
python_display_short_literal_types = True
python_use_unqualified_type_names = True


##############################################################################
#    Sphinx extension options                                                #
##############################################################################

# -- sphinx-bibtex (references) ----------------------------------------------
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "label"


# -- Docstring parsing (napoleon)---------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = True
napoleon_preprocess_types = False
napoleon_use_ivar = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = False


# -- Autodoc & autosummary ---------------------------------------------------
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": False,
    "special-members": "__call__",  #', __enter__, __iter__, __next__',
    "inherited-members": False,
    "show-inheritance": False,
    "member-order": "groupwise",
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"
autodoc_typehints_description_target = "documented"
autodoc_type_aliases = {
    "DataFrame": "pandas.DataFrame",
}
autoclass_content = "both"
autosummary_generate = glob.glob("*.rst")
autosummary_generate_overwrite = True


# -- Intersphinx -------------------------------------------------------------
# If you have trouble with proxies, etc., download these manually and place in
# the _local/ directory.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", (None, "_local/python-objects.inv")),
    "matplotlib": (
        "https://matplotlib.org/stable/",
        (None, "_local/matplotlib-objects.inv"),
    ),
    "numpy": ("https://numpy.org/doc/stable/", (None, "_local/numpy-objects.inv")),
    "pandas": ("https://pandas.pydata.org/docs/", (None, "_local/pandas-objects.inv")),
}


# -- Documenting extension modules with doxygen ------------------------------
if doxygen_installed:

    # -- Breathe options -----------------------------------------------------
    breathe_default_project = "libsansmic"
    breathe_separate_member_pages = False
    breathe_show_enumvalue_initializer = True
    breathe_domain_by_extension = {"hpp": "cpp"}
    breathe_projects_source = {
        "libsansmic": (
            "",
            glob.glob("../src/ext_modules/libsansmic/*"),
        )
    }
    breathe_default_members = ("members",)
    breathe_projects = {"libsansmic": "_build/doxyxml/xml"}
    breathe_show_include = False
    breathe_doxygen_aliases = {
        "rstref{1}": r"\verbatim embed:rst:inline :ref:`\1` \endverbatim"
    }
    breathe_order_parameters_first = True

    # -- Exhale options ------------------------------------------------------
    exhale_args = {
        "containmentFolder": "./apidocs",  # required
        "rootFileName": "index.rst",  # required
        "doxygenStripFromPath": "../src/ext_modules/libsansmic",  # required
        "rootFileTitle": "libsansmic",  # Heavily encouraged optional argument
        "createTreeView": False,
        "exhaleExecutesDoxygen": True,
        "exhaleDoxygenStdin": """INPUT = ../src/ext_modules/libsansmic
    EXCLUDE_SYMBOLS = "std"
    EXCLUDE_SYMBOLS += "PYBIND11_MODULE"
    SHOW_INCLUDE_FILES = "NO"
    SHOW_HEADER_FILES = "NO"
    HIDE_SCOPE_NAMES = "YES"
    EXTRACT_PRIVATE = "NO"
    EXTRACT_ALL = "NO"
    SHOW_USED_FILES = "NO"
    SHOW_FILES = "NO"
    SHOW_NAMESPACES = "NO"
    XML_PROGRAMLISTING = "NO"
    """,
        "contentsDirectives": False,
        "minifyTreeView": True,
        "fullToctreeMaxDepth": 2,
        "listingExclude": ["std", "PYBIND11_MODULE"],
        "unabridgedOrphanKinds": [],
    }


##############################################################################
#    Builder options                                                         #
##############################################################################

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
# html_sidebars = {"nomenclature": []}
html_theme_options = {
    "logo": {
        "image_light": "_static/logo-light.png",
        "image_dark": "_static/logo-dark.png",
    },
    "icon_links": [
        {
            "name": "GitHub",  # Label for this link
            "url": "https://github.com/sandialabs/sansmic",  # required URL where the link will redirect
            "type": "fontawesome",  # The type of image to be used
            "icon": "fa-brands fa-github",  # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
        },
        {
            "name": "Sandia National Laboratories",
            "url": "https://sandia.gov",  # required
            "type": "local",
            "icon": "_static/snl_logo.png",
        },
    ],
    "navigation_with_keys": False,
    "use_edit_page_button": False,
    "primary_sidebar_end": ["indices.html"],
    "show_toc_level": 2,
    "switcher": {
        "json_url": "https://sandialabs.github.io/sansmic/_static/switcher.json",
        "version_match": version,
    },
    # "secondary_sidebar_items": ["page-toc"], #["page-toc", "edit-this-page", "sourcelink"],
    "navbar_start": [
        "navbar-logo",
        "version-switcher",
    ],
    "navbar_end": [
        "theme-switcher",
        "navbar-icon-links",
    ],
    "analytics": {"google_analytics_id": "G-23TNKN36XM"},
}


# -- Options for LaTeX output -------------------------------------------------
latex_engine = "xelatex"
latex_use_xindy = False
latex_documents = [
    ("userman", "sansmic.tex", "User Guide for sansmic", "David Hart", "manual"),
]
latex_toplevel_sectioning = "chapter"
latex_domain_indices = False
latex_elements = {
    "papersize": r"letterpaper",
    "pointsize": r"10pt",
    "passoptionstopackages": r"""
\PassOptionsToPackage{svgnames}{xcolor}
""",
    "fontpkg": r"""
\usepackage{nimbussans}
\usepackage[nosf,nott]{kpfonts}
\usepackage[scaled=0.85]{sourcecodepro}
\usepackage[utf8x]{inputenc}
\usepackage[T1]{fontenc}
""",
    "preamble": r"""
\usepackage[titles]{tocloft}
\usepackage[
	range-units=single,
	per-mode=symbol,
	group-digits=decimal,
	group-minimum-digits=3,
	group-separator={,},
	range-phrase={\,--\,},
	quotient-mode=fraction,
	]{siunitx}\cftsetpnumwidth {1.25cm}
\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
\renewcommand\sphinxcrossref[1]{#1}
\renewcommand\sphinxtermref[1]{#1}
""",
    "sphinxsetup": "TitleColor=DarkGoldenrod",
    "fncychap": r"\usepackage[Bjornstrup]{fncychap}",
    "printindex": r"\footnotesize\raggedright\printindex",
}
