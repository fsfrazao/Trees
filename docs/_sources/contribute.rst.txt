==================
How to contribute
==================

*trees_ibm* is still in it's early stages of development and help is very much appreciated!

You can contribute by:

* Adding features to the codebase
* Expanding the testing suite
* Improving the documentation
* Reporting/fixing bugs
* Suggesting new features
* Giving your input on the biological meaning of the model assumptions/formulation




How to Report Issues
====================

To report an issue, please use the
`GitHub issues tracker <https://github.com/fsfrazao/Trees/issues>`_. When
reporting issues, please mention the following details:

* Which version you are using
* What was the source code that generated the problem (if applicable)
* Which platform are you running on
* How to reproduce the issue
* What was the result of the issue
* What were you expecting to get 

Reducing the source code that caused the issue to a bare minimum is always
very helpful.

Workflow for Pull Requests
==========================

In order to contribute, please fork off of the ``develop`` branch and make your
changes there. Your commit messages should detail *why* you made your change
in addition to *what* you did (unless it is a tiny change).

If you need to pull in any changes from ``develop`` after making your fork (for
example, to resolve potential merge conflicts), please avoid using ``git merge``
and instead, ``git rebase`` your branch.

Additionally, if you are writing a new feature, please ensure you write appropriate test cases and place them under ``src/tests/``.


Finally, please test your code and ensure that it runs locally before submitting a pull request.

Thank you for your help!

Running the tests
==========================

*trees-ibm* includes a set of tests. They are included in the /src/tests/  directory. They require the pytest package.

To run all tests go to the src directory and run: ::

  python -m pytest

You can also specify on module: ::

  python -m pytest test_world.py


