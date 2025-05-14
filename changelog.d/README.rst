Changelog
=========

This directory contains "news fragments" which are short files that contain a
small ReST-formatted text that will be added to the next ``CHANGELOG``.

Each fragment should concisely describe the changes made aimed for X-PSI users
and not developers.

Each file should be named like ``<pull_request_number>.<type>.rst``, where
``<pull_request_number>`` is a pull request number, and ``<type>`` is one of:

* ``summary``: A high-level overview of the main purpose or impact of the
  change.
* ``fixed``: Describes a bug fix or correction of faulty behavior.
* ``added``: Notes newly introduced features or functionality.
* ``changed``: Details modifications to existing behavior or interfaces.
* ``deprecated``: Marks features or APIs that are discouraged and will be
  removed in the future.
* ``removed``: Indicates features or components that have been taken out.
* ``attribution``: Credits contributors or acknowledges third-party work or
  inspiration.

For example: ``123.added.rst`` would have the content::

    The ``my_new_feature`` option is now available for ``my_favorite_function``.
    To use it, write ``np.my_favorite_function(..., my_new_feature=True)``.

Note the use of double-backticks for code.

If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

You can install ``towncrier`` and run ``towncrier --draft`` if you want to get a
preview of how your change will look in the final release notes.

.. note::

    This README was adapted from the numpy and pytest changelog readme under the
    terms of the MIT licence.