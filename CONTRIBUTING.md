# Contributing to sansmic

Thanks so much for your willingness to help out with sansmic's
development :grinning:  Here's everything you need to know.



## Contributor License Agreement

By contributing to this software project, you are agreeing to the following
terms and conditions for your contributions:

* You agree your contributions are submitted under the
  [BSD 3-Clause License](LICENSE).
* You represent you are authorized to make the contributions and grant the
  license.  If your employer has rights to intellectual property that includes
  your contributions, you represent that you have received permission to make
  contributions and grant the required license on behalf of that employer.



## Setting Up Your Environment

Before you begin hacking on sansmic, you'll need to do a few things:

1. **Fork the Repository:**  On [sansmic's main page][homepage], in
   the top-right corner, you'll see a button to [fork the repository][fork].
   Do so, and then clone your fork and set it up with
   ```bash
   git clone git@github.com:${USER}/sansmic
   cd sansmic
   git remote add upstream git@github.com:sandialabs/sansmic
   git fetch upstream
   ```
2. **Set Up your Environment:**  Install the requirements with
   ```bash
   pip install --user --editable .[dev,docs,test]
   ```
3. **Documentation:** To build the documentation, you will need to also
   obtain Doxygen as well as Sphinx. Sphinx will be installed with the
   pip command above, but you will need to [install Doxygen][doxygen].
   Afterwards, you can build a local copy of the documentation with
   ```bash
   cd docs
   make html
   ```

[homepage]: https://github.com/sandialabs/sansmic
[fork]: https://docs.github.com/en/get-started/quickstart/fork-a-repo
[doxygen]: https://www.doxygen.nl/

### Pre-Commit

We use [pre-commit][precommit] to ensure adherence to good software development
practices, enforce our style guide, etc.  To set yourself up with it, ensure
you have your development environment activated, and then run
```bash
cd /path/to/sansmic
pre-commit install
pre-commit install --hook-type commit-msg --hook-type pre-push --hook-type black
```

[precommit]: https://pre-commit.com/

The checks we perform are the following:
* Ensure no large files are added to the repository.
* Check for files that would conflict in case-sensitive filesystems.
* Ensure files don't contain merge conflict strings.
* Ensure files end with a single blank line.
* Ensure we only use Unix line endings.
* Trim trailing whitespace.
* Use black to format python code.

### VS Code

Our IDE of choice is [VS Code][vscode].  Download and install it, and then
follow the instructions below to get it set up.

[vscode]: https://code.visualstudio.com/

#### Getting Started

For getting up and running with VS Code, check out the following resources:

* [Visual Studio Code:  Introductory Videos][videos]
* [Visual Studio Code:  Getting Started][gettingstarted] documentation
* One-page reference of keyboard shortcuts for [Linux][linux], [MacOS][mac],
  and [Windows][win]

[videos]: https://code.visualstudio.com/docs/getstarted/introvideos
[gettingstarted]: https://code.visualstudio.com/docs
[linux]: https://code.visualstudio.com/shortcuts/keyboard-shortcuts-linux.pdf
[mac]: https://code.visualstudio.com/shortcuts/keyboard-shortcuts-macos.pdf
[win]: https://code.visualstudio.com/shortcuts/keyboard-shortcuts-windows.pdf

#### Extensions

VS Code is highly customizable through the use of extensions.  Click on the
building blocks icon on the left-hand side to open the **EXTENSIONS** panel to
search for and install them.  These are the ones we recommend:

**General**
* **Better Comments:**  Style comments in a more human-friendly manner.
* **Code Spell Checker:**  Catch spelling mistakes in code and comments.
* **Conventional Commits:**  Adhere to the [Conventional Commits][conventional]
  specification for commit messages.
* **Coverage Gutters:**  Display test coverage in the editor.
* **GitLens — Git supercharged:**  Integrate some of the powerful features of
  [GitKraken][kraken] into VS Code.
* **IntelliCode:**  AI-assisted development features.
* **Pre-Commit:**  Commands and helpers for executing pre-commit hooks.
* **Test Explorer UI:**  Extensible user interface for running your tests in VS
  Code.
* **Vim:**  For when you can't truly leave vi behind (and who would want to?).
* **vscode-icons:**  Icons for the file explorer.
* INSERT OTHERS
* AS NEEDED

[conventional]: https://www.conventionalcommits.org/en/v1.0.0/
[kraken]: https://www.gitkraken.com/

**Python-Specific**
* **Python:**  The basic extension for Python code in VS Code.
* **Pylance:**  Language server for syntax checking your python code.
* **isort:**  Sort imports in python code to be consistently ordered and prioritized.
* **Black formatter:**  Format python files without needing to wait for pre-commit.
* **Jupyter:**  Run and develop the example notebooks.

**C++-Specific**
* **C/C++:**  Basic C/C++ support.
* **C/C++ Extension Pack:**  Extra extensions that help with C++ development.

#### Settings

After installing the various extensions, you'll also want to customize your
**Settings**.  Some things can be set in the [Settings editor][settings]:

[settings]: https://code.visualstudio.com/docs/getstarted/settings#_settings-editor

* Text Editor
  * Font
    * **Font Family:**  Add "'Hack FC Ligatured', " to the beginning of the
      list.  You'll also want to install the [Ligatured Hack fonts][hack] on
      your system.
  * Files
    * **Auto Save:**  Set to "onFocusChange".
    * **Trim Trailing Whitespace:**  Check.
* Extensions
  * Conventional Commits
    * **Show Editor:**  Check.
    * **Silent Auto Commit:**  Check.
  * coverage-gutters
    * **Coverage-gutters: Show Line Coverage:**  Check.
    * **Coverage-gutters: Show Ruler Coverage:**  Check.
  * Git
    * **Allow Force Push:**  Check.
    * **Autofetch:**  Set to "true".
    * **Enable Smart Commit:**  Check.
    * **Fetch On Pull:**  Check.
    * **Prune On Fetch:**  Check.
    * **Rebase When Sync:**  Check.  This makes it such that you use `git pull
      --rebase` when pulling.
    * **Show Push Success Notification:**  Check.
    * **Terminal Git Editor:**  Check.
  * pre-commit-helper
    * **Run On Save:**  Select "all hooks".
  * Vim
    * **Vimrc:**  Enable:  Check.
  * Black
    * **Black:** Enable:  Check.
   

[hack]: https://github.com/gaplo917/Ligatured-Hack

Other items can be customized in the [`settings.json` file][settingsjson].
Consider adding the following snippets:

[settingsjson]: https://code.visualstudio.com/docs/getstarted/settings#_settingsjson

* To use ligatures to make multi-character symbols in code more readable, add
  ```json
  "editor.fontLigatures": true
  ```
* To add vertical rulers to the editor window, add
  ```json
  "editor.rulers": [
      72,
      79
  ]
  ```
* If you'd like the UI to preserve the scope that you're currently editing at
  the top of the file as you scroll through it, you can add
  ```json
  "editor.experimental.stickyScroll.enabled": true
  ```
* To turn on search highlighting when using the Vim extension, add
  ```json
  "vim.hlsearch": true
  ```

## Coding Standards

Many of our coding standards are applied and enforced using
[pre-commit][precommit] (see [above](#pre-commit)).  In addition to what's
handled automatically, we have the following.

### Documentation

All code should be appropriately documented using NumPy or Google style
docstrings (NumPy format is preferred) that can be parsed using Sphinx.
You should build the documentation locally to make sure that everything is
working correctly before creating a pull request.

Any new functionality that is added should also be added to the User Guide.
Documentation that goes into the User Guide must be appropriate, in English,
and free of typos and grammatical errors. If you are unsure about it, we
recommend copying your documentation into a word processor to verify that
it is correct. While we appreciate that spelling is different in different
countries, please use U.S. English spellings to keep documentation consistent
(sorry, but we can't modelling nor generating colour pictures here).

### Testing

All new functionality must come with unit tests. Existing tests must be run
and pass prior to a pull request being accepted. Do not change existing tests
without verifying with the maintainers that such action would be appropriate,
since doing so would indicate a major change in the codebase.

### Best Practices

Please create an issue, or cite an existing issue when you add code, as this
helps show that the codebase is moving forward with purpose.


## Creating Issues

Create [issues][issues] in GitHub for any work that needs to be done.

[issues]: https://github.com/sandialabs/sansmic/issues

### Markdown

[Markdown][markdown] is a lightweight markup language with plain text
formatting syntax.  GitHub uses a form of it for rendering issue and pull
request descriptions and comments, wiki pages, and any files in your
repositories with a `.md` extension (such as this one).  For more details on
what's possible with GitHub-flavored Markdown, [see the documentation][gfm].

[markdown]: https://en.wikipedia.org/wiki/Markdown
[gfm]: https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax

### Issue Templates

When creating an issue, select one of the available issue templates, and then
follow the instructions that appear in the template.

### Labels

The [sansmic > Labels][labels] page shows you all the labels we use,
along with their descriptions.  "Type" labels will be applied automatically
when an issue is created via the issue template.  "Stage" labels are used to
keep track of what work is being done at any given time.  Other labels can be
added as appropriate.

[labels]: https://github.com/sandialabs/sansmic/labels



## Working Issues

### Planning Work

Depending on the needs of the sansmic user base, a project maintainer
will schedule issues to be tackled in the near future by applying the
"Stage:  Soon" label.  An issue can then be grabbed by a member of the
development team when they have some time available to work on it.

> **Note:**  If you're working on something that hasn't been
> scheduled&mdash;e.g., you're tackling a quick bug fix, or are tinkering with
> a feature you'd like to see added&mdash;don't worry about the "Stage:  Soon"
> label; just work on it as per the rest of the guidelines below.

### When Work Begins

First switch the label from "Stage:  Soon" to "Stage:  Development".  Next make
sure your local `master` branch is up-to-date with
```bash
git checkout master
git pull --ff-only upstream master
```

> **Note:**  You should never be making commits on your `master` branch.  The
> `--ff-only` ensures you only update your local `master` branch if it can be
> fast-forwarded.

Once `master` is updated, you then create a feature branch off of it with
```bash
git checkout -b <branch-name>
```

The recommended branch naming convention is to use the issue number, followed
by a hyphen, followed by the issue title, all lowercase, omitting special
characters, and replacing spaces with hyphens.  For instance, if issue number
123 has "Implement Awesome New Feature" as the title, the corresponding branch
name would be `123-implement-awesome-new-feature`.

### As Work Continues

Do whatever work is necessary to address the issue you're tackling.  Break your
work into logical, working commits.  Use the **Conventional Commits** extension
for VS Code (or something similar) to ensure your commit messages adhere to the
[Conventional Commits specification][conventional].

Feel free to commit and push small chunks early and often and then use
interactive rebase to reorganize your commits before sharing.

> **Note:**  If you rebase a branch that's already been pushed to a remote,
> you'll wind up changing the history, which will require a force push.  That
> is permissible (even encouraged), but if you've had one or more reviewers or
> collaborators working with you on the branch, *get their buy-in first* before
> doing a force push.

### When Work is Complete

While working on your feature in your local `<branch-name>` branch, other
commits will likely make it into the upstream `master` branch.  There are a
variety of ways to merge these changes into your local feature branch.  One
possibility is
```bash
git checkout master
git pull --ff-only upstream master
git checkout <branch-name>
git rebase master
```

though there are others that are equally valid.  Once all is well, create a
pull request (see below).

### Closing Old Issues

If at any point you encounter an issue that will not be worked in the
foreseeable future, it is worthwhile to close the issue such that we can
maintain a reasonable backlog of upcoming work.  Do be sure to include in the
comments some explanation as to why the issue won't be addressed.



## Pull Requests

The only way changes get into `master` is through pull requests.  When you've
completed work on an issue, push your branch to your fork with `git push -u
origin <branch-name>`, and then create a pull request.  Apply the same "Type"
label as the issue, and then return to the issue and swap "Stage:  Development"
for "Stage:  Review".

### Reviews

A project maintainer will review your pull request.  Work with them to get your
changes into an acceptable state.

### Drafts

You may wish to have your changes reviewed by colleagues before they are ready
to be merged into `master`.  To do so, create a [draft pull request][drafts].
GitHub will not allow you to merge a draft request.  When it's ready for
review, you can [mark it as ready][ready].

[drafts]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests#draft-pull-requests
[ready]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/changing-the-stage-of-a-pull-request#marking-a-pull-request-as-ready-for-review

### Continuous Integration

GitHub Actions is GitHub's continuous integration and delivery (CI/CD)
mechanism.  The configurations for our workflows can be found in the `*.yml`
files in [the `.github/workflows` directory](.github/workflows) in the
repository root:
* INSERT

In addition to the GitHub Actions workflows, we also integrate with the
following services:
* [Codecov][codecov]:  We calculate and publish code coverage details for all
  PRs and merged into `master`, and the service comments on PRs indicating any
  changes in coverage compared to the base branch.
* [pre-commit.ci][precommitci]: We ensure all the `pre-commit` checks we
  encourage developers to run locally also run in CI.
* INSERT OTHERS

[codecov]: https://about.codecov.io/
[precommitci]: https://pre-commit.ci/

### Merging

When the review is finished and changes are ready to be merged into `master`:
1. Rebase your feature branch on top of the latest `master`.
2. Clean up your branch with an interactive rebase, squashing down to the
   smallest number of commits that makes sense.  If there are successive
   distinct changes in a pull request, it's fine for those to be preserved in
   separate commits.
3. Notify your reviewers that the request is ready to merge.
4. Wait for them to merge the request.
5. Ensure the post-merge CI pipeline also succeeds.



## Conventional Comments

When commenting on issues and pull requests, we endeavor to adhere to the
[Conventional Comments specification][conventionalcomments] wherever
applicable.  Comments should have the form
```
<label> [decorations]: <subject>

[discussion]
```
where:
* **label:**  This is a single label that signifies what kind of comment is
  being left.
* **subject:**  This is the main message of the comment.
* **decorations (optional):**  These are extra decorating labels for the
  comment.  They are surrounded by parentheses and comma-separated.
* **discussion (optional):**  This contains supporting statements, context,
  reasoning, and anything else to help communicate the “why” and “next steps”
  for resolving the comment.

The labels we use are the following:

* **chore:**  Chores are simple tasks that must be done before the subject can
  be “officially” accepted.  Usually, these comments reference some common
  process.  Try to leave a link to the process description so that the reader
  knows how to resolve the chore.
* **issue:**  Issues highlight specific problems with the subject under review.
  These problems can be user-facing or behind the scenes.  It is strongly
  recommended to pair this comment with a suggestion.  If you are not sure if a
  problem exists or not, consider leaving a **question**.
* **note:**  Notes are always non-blocking and simply highlight something the
  reader should take note of.
* **polish:**  Polish comments are like a suggestion, where there is nothing
  necessarily wrong with the relevant content, there's just some ways to
  immediately improve the quality.
* **praise:**  Praises highlight something positive.  Try to leave at least one
  of these comments per review.  Do not leave false praise (which can actually
  be damaging).  Do look for something to sincerely praise.
* **quibble:**  Quibbles are trivial, preference-based requests.  These should
  be non-blocking by nature.
* **suggestion:**  Suggestions propose improvements to the current subject.
  It's important to be explicit and clear on what is being suggested and why it
  is an improvement.  Consider using patches and the blocking or non-blocking
  decorations to further communicate your intent.
* **thought:**  Thoughts represent an idea that popped up from reviewing.
  These comments are non-blocking by nature, but they are extremely valuable
  and can lead to more focused initiatives and mentoring opportunities.
* **todo:**  TODO's are small, trivial, but necessary changes.  Distinguishing
  **todo** comments from **issues** or **suggestions** helps direct the
  reader's attention to comments requiring more involvement.
* **typo:**  Typo comments are like **todo**, where the main issue is a
  misspelling.
* **question:**  Questions are appropriate if you have a potential concern but
  are not quite sure if it's relevant or not.  Asking the author for
  clarification or investigation can lead to a quick resolution.

The decorations we use are the following:

* **if-minor:**  This decoration gives some freedom to the author that they
  should resolve the comment only if the changes ends up being minor or
  trivial.
* **non-blocking:**  A comment with this decoration should not prevent the
  subject under review from being accepted.  Without this decoration, comments
  are assumed to be blocking.

[conventionalcomments]: https://conventionalcomments.org/
