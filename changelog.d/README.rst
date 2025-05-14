Changelog
=========

This directory contains "news fragments" which are short files that contain a small ReST-formatted text that will be added to the next ``CHANGELOG``.

Each fragment should concisely describe the changes made aimed for X-PSI users and not developers.

Each file should be named like ``<pull_request_number>.<type>.rst``, where
``<pull_request_number>`` is a pull request number, and ``<type>`` is one of:

* ``summary``: A change which is not backwards compatible and requires users to change code.
* ``fixed``: New user facing features and any new behavior.
* ``added``: Fixes a reported bug.
* ``changed``: Documentation improvement, like rewording an entire session or adding missing docs.
* ``deprecated``: Feature deprecation.
* ``removed``: Feature removal.
* ``attribution``: Fixes a small typo or internal change that might be noteworthy.

For example: ``123.added.rst`` would have the content::

    The ``my_new_feature`` option is now available for ``my_favorite_function``.
    To use it, write ``np.my_favorite_function(..., my_new_feature=True)``.

Note the use of double-backticks for code.

If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

You can install ``towncrier`` and run ``towncrier --draft`` if you want to get a preview of how your change will look in the final release
notes.

.. note::

    This README was adapted from the numpy and pytest changelog readme under the terms of the MIT licence.