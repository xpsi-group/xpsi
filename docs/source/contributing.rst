.. _contributing:

Contributing
------------

We outline below the workflow we use for developing X-PSI.

Contact
~~~~~~~

You are welcome to contact members of the X-PSI team directly if you have
any questions about the software or its use (for the current active list see
the :ref:`acknowledgements` page).
To get in touch with us regarding bugs and issues the easiest way is via the
GitHub `Issues <https://github.com/xpsi-group/xpsi/issues/>`_ page.

Contributing
~~~~~~~~~~~~

To contribute to this project you can do so whilst being recognised as either
a community member or a team member. If you contribute to feature development
publicly on GitHub, you may well be invited to be recognised as a team member in
the :ref:`acknowledgements`.

Past developers are recognised as team members unless they no longer wish to be.
Aside from acknowledging your contribution, this is useful for feature
maintenance and as a contact point for questions, if permission is granted.
Past team members will be acknowledged as such.

If you use generative AI tools (such as large language models or AI coding
assistants) while preparing a contribution, please read our :ref:`ai_policy`.
In short: you are responsible for everything you submit, every pull request
must include an AI disclosure (stating the tools, models and versions used, if
any), and we do not accept contributions opened by automated agents.

If you wish to base a major research project on the
development and/or application of X-PSI, we invite you (and your
group where applicable) to discuss joining the X-PSI team in order to
collaborate. If you find this idea interesting then please contact one of the core team members
(contact details are on the :ref:`acknowledgements` page).


.. _workflow:

Git workflow
~~~~~~~~~~~~

Advisory workflow for contributing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

X-PSI is an open-source project and we aim to ensure that all X-PSI versions required
for reproducibility of publications have release tags on GitHub.

End-users and community contributors can interact with the project freely on
GitHub. If you wish to contribute code, either in the form of a patch, a tweak,
or a feature, but you are uncertain about the workflow, we provide the advisory
workflow below. It follows GitHub's standard
`fork and pull request
<https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-open-source>`_
model. Throughout, ``origin`` refers to your fork and ``upstream`` to the
central ``xpsi-group/xpsi`` repository.

**Set up (once)**

* Create your own fork of the ``xpsi`` repository by clicking **Fork** on the
  `X-PSI GitHub page <https://github.com/xpsi-group/xpsi>`_.

* Clone your fork to your local computer. Git names this remote ``origin``:

.. code-block:: bash

    git clone https://github.com/<username>/xpsi.git <path/to/xpsi>
    cd <path/to/xpsi>

* Add the central repository as a second remote called ``upstream``, and check
  that both remotes are set:

.. code-block:: bash

    git remote add upstream https://github.com/xpsi-group/xpsi.git
    git remote -v

* Install X-PSI from your clone as described in :ref:`install`. X-PSI contains
  compiled extensions, so re-install (``pip install .``) after changing them.

.. note::

    Older versions of this page told you to clone ``xpsi-group/xpsi`` first (as
    ``origin``) and to add your fork later as ``fork``. If your clone is set up
    that way you can keep using it, replacing ``upstream`` with ``origin`` and
    ``origin`` with ``fork`` in the commands below, or rename the remotes to
    match this page:

    .. code-block:: bash

        git remote rename origin upstream
        git remote rename fork origin

**Keep your main branch up to date**

* Use your ``main`` branch only as a copy of ``upstream/main``, and do not
  commit to it directly. That way, new upstream commits can always be applied to
  ``main`` with a fast-forward merge:

.. code-block:: bash

    git checkout main
    git fetch upstream
    git merge --ff-only upstream/main
    git push origin main

* Alternatively, click **Sync fork** on your fork's GitHub page (or run
  ``gh repo sync`` with the GitHub CLI) and then run ``git pull`` on your
  local ``main``. See GitHub's guide to
  `syncing a fork
  <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork>`_.

**Make your changes**

* Starting from an up-to-date ``main``, create a branch for your work, which
  we'll assume is some patch, but could be a feature or otherwise:

.. code-block:: bash

    git checkout -b patch/fix_something

* Commit your work to ``patch/fix_something``:

.. code-block:: bash

    git add <changed files>
    git commit -m 'Fix something'

* Before opening a pull request, run the test suite and the checks that our
  continuous integration (CI) runs on every pull request:

.. code-block:: bash

    pytest
    flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics

* Add tests and documentation for new or changed functionality where
  appropriate.

**Integrate upstream changes**

* Check whether there are any new upstream commits, which would mean that
  ``patch/fix_something`` and ``upstream/main`` have diverged:

.. code-block:: bash

    git fetch upstream

* If the branches have diverged, you can either rebase ``patch/fix_something``
  on ``upstream/main`` or merge ``upstream/main`` into
  ``patch/fix_something``, in either case resolving any conflicts:

.. code-block:: bash

    git rebase upstream/main <or> git merge upstream/main

* Once you have pushed ``patch/fix_something`` to your fork, and especially
  once you have opened a pull request from it, only merge (do not rebase), so
  that the branch history is preserved and you do not need to force push. (The
  one exception is fixing commits flagged by the AI co-author check described
  in our :ref:`ai_policy`, which needs a force push with
  ``--force-with-lease``, only on your own pull request branch and never on
  ``main``; the check explains what to do.)

* After integrating upstream changes, you might decide to continue working on
  your branch. Alternatively, you could work on another branch for a different
  patch or feature. In either case you should repeat the process of integrating
  upstream changes as appropriate, and as a requirement if preparing for a pull
  request (although there could be exceptional corner cases wherein an X-PSI
  team member assists with this merging process).

**Open a pull request**

* Push ``patch/fix_something`` to your fork, creating a remote branch
  ``origin/patch/fix_something`` that ``patch/fix_something`` tracks:

.. code-block:: bash

    git push -u origin patch/fix_something

* On GitHub, open a pull request from ``patch/fix_something`` on your fork to
  ``main`` on ``xpsi-group/xpsi`` (GitHub usually offers a **Compare & pull
  request** button after you push). In the description:

  - explain what the pull request changes and why;
  - reference any relevant issues, writing for example ``Closes #123`` for an
    issue that the pull request fixes, so that GitHub closes it automatically
    when the pull request is merged;
  - include an AI disclosure, as described in our :ref:`ai_policy`. This is
    required for every pull request, including those where no AI tools were
    used.

* Leave **Allow edits by maintainers** ticked, so that X-PSI team members can
  push small fixes to your branch if needed.

* If your work is not yet ready for review but you would like feedback, open
  it as a **draft pull request**.

* Add changelog fragments, which are short files in the ``changelog.d/``
  directory named ``<pull request number>.<type>.rst``, for example
  ``123.fixed.rst``. Since they use the pull request number, add them after
  opening the pull request and push them to the same branch. Include an
  ``attribution`` fragment with your name. The available types and formatting
  are explained in `changelog.d/README.rst
  <https://github.com/xpsi-group/xpsi/blob/main/changelog.d/README.rst>`_.
  The CI checks every pull request for these fragments; for very small changes
  a team member can add the ``no-changelog-needed`` label instead.

* You can update the pull request by pushing additional commits from
  ``patch/fix_something`` to your fork, which will update the pull request
  automatically:

.. code-block:: bash

    git push

Pull Request Handling and Development Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Pull request review and outcome:** The pull request will then be reviewed
  and discussed. The proposed changes will either be merged or a merge will be
  pending because we request additional commits on the pull-request topic
  branch. Your pull request may be declined in some instances because the work
  reproduces development work that has already been performed but not published;
  your pull request may also be ultimately declined if it contains changes or
  implementations that we do not support or cannot maintain, and which cannot
  for some reason be separated from changes we do support and could maintain.
  Your intellectual contribution to the project will be gratefully acknowledged
  in the :ref:`acknowledgements` and/or in the project :ref:`history` if this
  interaction leads to some form of merged development/implementation by another
  community member, even if your pull request is ultimately declined.

* **Acknowledging co-authors:** If you co-authored a pull request with one or
  more collaborators, you can acknowledge them using the GitHub pull-request GUI
  as you would for a single commit. When a pull request is accepted, it is
  typically going to be via a merge-squash unless the history is clean or work
  will continue to be committed on the topic branch after the merge (where
  applicable). In this case it is the responsibility of the X-PSI team member
  executing the merge to replicate the list of co-authors from the original pull
  request in the squash message. Co-authors are people: AI tools should not be
  listed as co-authors (see :ref:`ai_policy`).

* **Handling merge conflicts:** When a pull request is merged, conflicts will
  either need to be resolved locally by you as suggested above, ending in a pull
  request update, or by an X-PSI team member locally and then merged with or
  without a pull request.

* **Developing new features or patches:** If you are ready to start development
  on a distinct patch or feature that is not conditional on your open pull
  requests being merged, then you can apply the workflow above by branching
  (again) off of an up-to-date ``main``. If your work *is* conditional on your
  open pull requests, you are free to continue your development by committing to
  the relevant topic branch (or according to some other branching scheme).
  However, there is a risk that more work will be needed if the open pull
  request is not merged into the central repository; or if only a subset of
  proposed changes are merged; or conflict resolution does not favour all of the
  changes you proposed. Of course, such work may nevertheless remain useful in
  your own applications even if it is never all merged into the central
  repository. If the pull request is merged after your continuation, and the
  plan is to submit a future pull request, you will have to merge in the
  ``upstream/main`` branch before opening another pull request so that the merge
  conflicts that were already resolved are not raised again.

* **Working with non-main branches:** The above workflow also applies to remote
  branches other than ``main`` that might exist in the ``xpsi`` repository that
  you wish to contribute to, but this should be a less common pattern.

If you want to contribute a feature, you are welcome to communicate with us
either on GitHub via issues and pull-requests, or by contacting a team member
directly.
