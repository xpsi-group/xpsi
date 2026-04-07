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
github `Issues <https://github.com/xpsi-group/xpsi/issues/>`_ page. 

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

If you wish to base a major research project on the
development and/or application of X-PSI, we invite you to discuss you (and your
group where applicable) joining the X-PSI team in order to collaborate. If you
find this idea interesting then please contact Anna Watts (contact details
are on the :ref:`acknowledgements` page).


.. _workflow:

Git workflow
~~~~~~~~~~~~

Advisory workflow for contributing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

X-PSI is an open-source project and we aspire that all X-PSI versions required
for reproducibility of publications have release tags on GitHub.

End-users and community contributors can interact with the project freely on
GitHub. If you wish to contribute code, either in the form of a patch, a tweak,
or a feature, but you are uncertain about the workflow, we now provide an
advisory workflow.

* Clone the ``xpsi`` repository to your local computer:

.. code-block:: bash

    git clone https://github.com/xpsi-group/xpsi.git <path/to/xpsi>

* You will be on the ``main`` branch by default; this branch by default tracks
  the official (or central) ``origin/main`` branch. Moreover, ``main``
  is merely to act as an image of ``origin/main``, meaning that official
  upstream commits can always be applied to ``main`` via a fast-foward merge.

* Checkout a local branch for your work, which we'll assume is some patch, but
  could be a feature or otherwise:

.. code-block:: bash

    git checkout -b patch/fix_something

* Commit your work to ``patch/fix_something``:

.. code-block:: bash

    git commit -m 'patch something'

* Check to see if any there are any new upstream commits, which would mean
  that ``patch/fix_something`` and ``origin/main`` have diverged:

.. code-block:: bash

    git fetch origin main:main

* If the branches have diverged, you can either rebase ``patch/fix_something``
  on ``main`` or merge ``main`` into ``patch/fix_something``, in either case
  resolving any conflicts:

.. code-block:: bash

    git rebase main <or> git merge main

* Note that if you have already pushed ``patch/fix_something`` to a remote
  repository you own (such as a fork of ``xpsi``; see below), and especially if
  this is accessible by others (e.g., via  submitted pull request), you should
  only consider merging ``main`` into ``patch/fix_something`` in order to
  preserve the branch history.

* After integrating upstream changes, you might decide to continue working on
  your branch. Alternatively, you could work on another branch for a different
  patch or feature. In either case you should repeat the process of integrating
  upstream changes as appropriate, and as a requirement if preparing for a pull
  request (although there could be exceptional corner cases wherein an X-PSI
  team member assists with this merging process).

* Once you are ready to contribute your work to the ``xpsi`` repository,
  generally meaning that you have integrated any upsteam changes from ``xpsi``,
  you need a fork of the ``xpsi`` repository on the same hosting platform
  (GitHub). You can create a fork using the GitHub GUI.

* With the address of your ``fork`` you can add it as a remote to your local
  repository:

.. code-block:: bash

    git remote add fork https://github.com/<username>/xpsi.git

* Now push ``patch/fix_something`` to ``fork``, creating a remote branch
  ``fork/patch/fix_something`` that ``patch/fix_something`` tracks:

.. code-block:: bash

    git push -u fork

* Now you can submit a pull request, using the GitHub GUI, from
  ``fork/patch/fix_something`` to ``xpsi/main``. Please reference any open
  issues that are to be closed or are relevant to the proposed changes.

* You can update the pull-request topic branch by pushing additional commits
  from ``patch/fix_something`` to ``fork/patch/fix_something``, which will
  update the pull request automatically:

.. code-block:: bash

    git push

Pull Request Handling and Development Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Pull request review and outcome:** The pull request will then be reviewed
  and discussed. The proposed changes will either be  merged or a merge will be
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
  will continue to be commited on the topic branch after the merge (where
  applicable). In this case it is the responsiblity of the X-PSI team member
  executing the merge to replicate the list of co-authors from the original pull
  request in the squash message.

* **Handling merge conflicts:** When a pull request is merged, conflicts will
  either need to be resolved locally by you as suggested above, ending in a pull
  request update, or by an X-PSI team member locally and then merged with or
  without a pull request.

* **Developing new features or patches:** If you are ready to start development
  on a distinct patch or feature that is not conditional on your open pull
  requests being merged, then you can apply the workflow above by branching
  (again) off of an up-to-date ``main``. If your work *is* conditional on your
  open pull requests, you are free to continue your development by commiting to
  the relevant topic branch (or according to some other branching scheme).
  However, there is a risk that more work will be needed if the open pull
  request is not merged into the central repository; or if only a subset of
  proposed changes are merged; or conflict resolution does not favour all of the
  changes you proposed. Of course, such work may nevertheless remain useful in
  your own applications even if it is never all merged into the central
  repository. If the pull request is merged after your continuation, and the
  plan is submit a future pull request, you will have to merge in the
  ``xpsi/main`` branch before opening another pull request so that the merge
  conflicts that were already resolved are not raised again.

* **Working with non-main branches:** The above workflow also applies to remote
  branches other than ``main`` that might exist in the ``xpsi`` repository that
  you wish to contribute to, but this should be a less common pattern.

If you want to contribute a feature, you are welcome to communicate with us
either on GitHub via issues and pull-requests, or by contacting a team member
directly. 
